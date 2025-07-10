from __future__ import annotations

from typing import Any, TypedDict
from pathlib import Path
from textwrap import dedent
import json

import pandas as pd
import scanpy as sc
from anndata import AnnData
from matplotlib import pyplot as plt

from latch_curate.utils import _fig_to_base64, write_anndata, write_html_report, df_to_str
from latch_curate.constants import latch_curate_constants as lcc
from latch_curate.utils import _df_to_html
from latch_curate.llm_utils import prompt_model
from latch_curate.tinyrequests import post
from latch_curate.config import user_config

statistics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]

class StatisticInterval(TypedDict):
    statistic_name: str
    interval:  tuple[float, float]
    reasoning: str

def qcol(metric: str, q: float) -> str:
    return f"{metric}_q{int(q*1000):03d}"

def _violin_plot(adata: AnnData, groupby: str | None = None):
    if adata.n_obs == 0 or (groupby and adata.obs[groupby].nunique() == 0):
        return plt.figure()
    if groupby and pd.api.types.is_categorical_dtype(adata.obs[groupby]):
        adata.obs[groupby] = adata.obs[groupby].cat.remove_unused_categories()
    grid = sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        groupby=groupby,
        rotation=45 if groupby else 0,
        show=False,
    )
    if isinstance(grid, list):
        return grid[0].figure
    return grid.fig


def compute_quantiles(adata: AnnData) -> pd.DataFrame:
    tail_quants = {"low": [0.01, 0.05], "high": [0.90, 0.95, 0.99, 0.995]}
    metrics = {
        "n_genes_by_counts": ["low", "high"],
        "total_counts":      ["low", "high"],
        "pct_counts_mt":     ["high"],
    }

    return (
        adata.obs.groupby("latch_sample_id").apply(
            lambda df: pd.Series(
                {
                    f"{m}_q{int(q*1000):03d}": df[m].quantile(q)
                    for m, tails in metrics.items()
                    for t in tails
                    for q in tail_quants[t]
                }
            )
        )
    )
def build_adaptive_threshold_prompt(
    stat_names: list[str],
    quantile_table: pd.DataFrame,
    sample_list: list[str],
) -> str:
    sample_list = sorted(sample_list)

    example_output = {
        sid: {
            "n_genes_by_counts": {"interval": [20, 6000], "reasoning": "keep central 98 % of genes"},
            "total_counts":      {"interval": [None, 2_000], "reasoning": "≤99.5 th percentile of UMIs"},
            "pct_counts_mt":     {"interval": [None, 15],    "reasoning": "cells with high MT are stressed"},
        }
        for sid in sample_list
    }

    return dedent(f"""
        You are an expert single-cell QC assistant.  Your task is to define
        per-sample filtering *intervals* for the following statistics:

          • n_genes_by_counts (total detected genes per cell)  
            – use BOTH lower and upper bounds  
          • total_counts (total UMI counts per cell)  
            – usually only an UPPER bound; set lower bound to null if unused  
          • pct_counts_mt (percentage of mitochondrial UMIs)  
            – usually only an UPPER bound; set lower bound to null if unused

        <quantile_table>
        {quantile_table.to_csv()}
        </quantile_table>

        Rules
          1. Derive sensible bounds from the quantiles of each sample
             (e.g. q01/q99, q05/q995, etc.).  
          2. Return a JSON object whose keys are the sample IDs and whose
             values are dictionaries exactly like in the example below.  
          3. Bounds that do not apply **MUST be null**.  
          4. Provide concise scientific reasoning (<40 words).  
          5. Output *only* the raw JSON – no markdown, no commentary.

        Example of the required structure (values are illustrative):

        {json.dumps(example_output, indent=2)}
    """)

def apply_fixed_thresholds(adata: AnnData, min_genes: int, max_counts: int, max_pct_mt: int):
    sc.pp.filter_cells(adata, min_genes=min_genes)
    if max_counts is not None:
        sc.pp.filter_cells(adata, max_counts=max_counts)
    adata = adata[adata.obs.pct_counts_mt < max_pct_mt]
    return adata

def build_adaptive_mask(
    adata: AnnData,
    interval_data: dict[str, list[StatisticInterval]],
) -> pd.Series:

    global_mask = pd.Series(False, index=adata.obs.index, dtype=bool)

    for sample_id, stat_list in interval_data.items():
        sample_mask = adata.obs["latch_sample_id"] == sample_id

        for stat in stat_list:
            stat_name = stat["statistic_name"]

            lower, upper = stat["interval"]

            if lower is not None:
                sample_mask &= adata.obs[stat_name] >= lower
            if upper is not None:
                sample_mask &= adata.obs[stat_name] <= upper

        global_mask |= sample_mask

    return global_mask


def build_qc_report_html(stages: list[dict[str, Any]], title: str = "QC & Filtering Report") -> str:
    rows_html: list[str] = []
    for s in stages:
        rows_html.append(
            f"<h2>{s['name']}</h2>"
            f"<p><strong>Total cells:</strong> {s['total']}</p>"
            f"<h3>Per-sample counts</h3>\n<table><tr><th>Sample ID</th><th>Count</th></tr>"
            + "".join(f"<tr><td>{sid}</td><td>{c}</td></tr>" for sid, c in s["counts"].items())
            + "</table>"
            + "<h3>Quantiles</h3>" + s["qtable"]
            + "<h3>Distributions</h3>"
            + f"<img src='data:image/png;base64,{s['global_fig']}' alt='{s['name']} global'/>"
            + f"<img src='data:image/png;base64,{s['sample_fig']}' alt='{s['name']} by-sample'/>"
            + (
                "<h3>Adaptive intervals &amp; reasoning</h3>" + s["interval_html"]
                if s["interval_html"]
                else ""
            )
        )

    body = "".join(rows_html)
    return dedent(
        f"""<!DOCTYPE html>
        <html><head><meta charset='utf-8'><title>{title}</title>
        <style>
          body {{ font-family: sans-serif; max-width: 900px; margin:auto; }}
          h2 {{ border-bottom: 1px solid #ccc; margin-top:2em; }}
          table {{ border-collapse: collapse; width:100%; margin-bottom:1em; }}
          th,td {{ border:1px solid #ddd; padding:4px 8px; text-align:center; }}
          img {{ max-width:100%; height:auto; margin-bottom:1em; }}
          .quant {{ font-size: 0.8em; }}
        </style></head><body>
          <h1>{title}</h1>
          {body}
        </body></html>"""
    )

def _intervals_to_html(intervals: dict[str, list[StatisticInterval]]) -> str:
    if not intervals:
        return ""
    rows = []
    for sample, stats in intervals.items():
        rows.append(f"<h4>Sample {sample}</h4>")
        rows.append(
            "<table>"
            "<tr><th>Statistic</th><th>Lower</th><th>Upper</th><th>Reasoning</th></tr>"
            + "".join(
                f"<tr><td>{s['statistic_name']}</td>"
                f"<td>{s['interval'][0]}</td>"
                f"<td>{s['interval'][1]}</td>"
                f"<td class='reason'>{s['reasoning']}</td></tr>"
                for s in stats
            )
            + "</table>"
        )
    return "".join(rows)

def _stage_snapshot(
    name: str,
    _adata: AnnData,
    *,
    intervals: dict[str, list[StatisticInterval]] | None = None,
) -> dict[str, Any]:
    return {
        "name": name,
        "total": _adata.n_obs,
        "counts": _adata.obs["latch_sample_id"].value_counts().to_dict(),
        "global_fig": _fig_to_base64(_violin_plot(_adata)),
        "sample_fig": _fig_to_base64(_violin_plot(_adata, groupby="latch_sample_id")),
        "qtable": _df_to_html(compute_quantiles(_adata)),
        "interval_html": _intervals_to_html(intervals or {}),
    }


def qc_and_filter(
    adata: AnnData,
    study_metadata_path: Path,
    paper_text_path: Path,
    workdir: Path,
    use_params: bool
):
    params_file = (workdir / lcc.qc_params_name)
    if use_params:
        assert params_file.exists(), "Missing parameters file"
        with open(params_file) as f:
            params = json.load(f)

    assert 'latch_sample_id' in adata.obs.columns
    adata.obs['latch_sample_id'] = adata.obs['latch_sample_id'].astype('category')

    workdir.mkdir(exist_ok=True)

    print("Calculating qc metrics: mt‑genes, n_genes_by_counts, total_counts, pct_counts_mt")

    adata.var["mt"] = adata.var["gene_symbols"].str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=False, inplace=True)

    study_metadata = study_metadata_path.read_text()
    paper_text = paper_text_path.read_text()
    quantile_table = compute_quantiles(adata)

    if use_params:
        fixed = params['fixed']
        min_genes = fixed['min_genes']
        max_counts = fixed['max_counts']
        max_pct_mt = fixed['max_pct_mt']
        print("Using existing fixed thresholds from params file.")
    else:
        print("Requesting fixed thresholds from model")
        resp = post(
            f"{lcc.nucleus_url}/{lcc.get_fixed_qc_thresholds_endpoint}",
            {
                "study_metadata": study_metadata,
                "paper_text": paper_text,
                "quantile_table": df_to_str(quantile_table),
                "metadata": json.dumps({"step": "qc", "project":
                                        workdir.resolve().parent.name}),
                "session_id": -1
            },
            headers = {"Authorization": f"Latch-SDK-Token {user_config.token}"}
        )
        try:
            data = resp.json()['data']
            min_genes = data['min_genes']
            max_counts = data['max_counts']
            max_pct_mt = data['max_pct_mt']
        except KeyError:
            raise ValueError(f'Malformed response data :{data}')

    stages: list[dict[str, Any]] = []
    stages.append(_stage_snapshot("Before filtering", adata))

    print("Applying fixed thresholds")

    adata = apply_fixed_thresholds(adata, min_genes, max_counts, max_pct_mt)

    stages.append(_stage_snapshot("After fixed thresholds", adata))

    print("Applying adaptive thresholds")

    interval_data: dict[str, list[StatisticInterval]] = {}

    if use_params:
        interval_data = params['adaptive']
        print("Using existing adaptive thresholds from params file")
    else:
        adaptive_threshold_prompt = build_adaptive_threshold_prompt(statistics, quantile_table, set(adata.obs['latch_sample_id']))
        while True:
            print("Requesting adaptive thresholds from language model")
            message_resp_json, _ = prompt_model([{"role": "user", "content": adaptive_threshold_prompt}])
            try:
                data = json.loads(message_resp_json)
                for sample_name, x in data.items():

                    for s in statistics:
                        assert s in x

                    interval_data[sample_name] = []
                    for k, v in x.items():
                        assert "interval" in v
                        assert "reasoning" in v
                        stat_data = StatisticInterval(statistic_name=k, interval=v["interval"], reasoning=v["reasoning"])
                        interval_data[sample_name].append(stat_data)
                break
            except Exception as e:
                print(e.with_traceback())
                print(f"Invalid model response: {e}. Trying again")
                continue

    adaptive_mask = build_adaptive_mask(adata, interval_data)

    print(f"Adaptive filter retains {adaptive_mask.sum()} / {adaptive_mask.size} cells")
    adata = adata[adaptive_mask, :]

    stages.append(_stage_snapshot("After adaptive thresholds", adata, intervals=interval_data))

    write_anndata(adata, workdir, lcc.qc_adata_name)

    html = build_qc_report_html(stages)
    write_html_report(html, workdir, lcc.qc_report_name)

    if not use_params:
        param_data = {
                "fixed": {"min_genes": min_genes, "max_counts": max_counts, "max_pct_mt": max_pct_mt},
                "adaptive": interval_data
        }
        with open(workdir / lcc.qc_params_name, "w") as f:
            json.dump(param_data, f)
        print(f"QC parameters written to {workdir / lcc.qc_params_name}")
    print(f"QC pipeline finished successfully – final cell count: {adata.n_obs}")

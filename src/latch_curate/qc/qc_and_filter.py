from __future__ import annotations

from typing import Any
from pathlib import Path
from textwrap import dedent
import json

import pandas as pd
import scanpy as sc
from anndata import AnnData
from matplotlib import pyplot as plt

from latch_curate.tinyrequests import post
from latch_curate.utils import _fig_to_base64, write_anndata, write_html_report
from latch_curate.constants import latch_curate_constants as lcc
from latch_curate.config import user_config
from latch_curate.utils import _df_to_html, df_to_str

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

def build_qc_report_html(stages: list[dict[str, Any]], title: str = "QC & Filtering Report") -> str:
    rows_html: list[str] = []
    for s in stages:
        rows_html.append(
            f"<h2>{s['name']}</h2>"
            f"<p><strong>Total cells:</strong> {s['total']}</p>"
            f"<h3>Per‑sample counts</h3>\n<table><tr><th>Sample ID</th><th>Count</th></tr>"
            + "".join(f"<tr><td>{sid}</td><td>{c}</td></tr>" for sid, c in s["counts"].items())
            + "</table>"
            + "<h3>Quantiles</h3>" + s["qtable"]
            + "<h3>Distributions</h3>"
            + f"<img src='data:image/png;base64,{s['global_fig']}' alt='{s['name']} global'/>"
            + f"<img src='data:image/png;base64,{s['sample_fig']}' alt='{s['name']} by‑sample'/>"
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

def _stage_snapshot(name: str, _adata: AnnData) -> dict[str, Any]:
    return {
        "name": name,
        "total": _adata.n_obs,
        "counts": _adata.obs["latch_sample_id"].value_counts().to_dict(),
        "global_fig": _fig_to_base64(_violin_plot(_adata)),
        "sample_fig": _fig_to_base64(_violin_plot(_adata, groupby="latch_sample_id")),
        "qtable": _df_to_html(compute_quantiles(_adata)),
    }


def qc_and_filter(
    adata: AnnData,
    study_metadata_path: Path,
    paper_text_path: Path,
    workdir: Path
):
    assert 'latch_sample_id' in adata.obs.columns
    adata.obs['latch_sample_id'] = adata.obs['latch_sample_id'].astype('category')

    workdir.mkdir(exist_ok=True)

    print("Calculating qc metrics: mt‑genes, n_genes_by_counts, total_counts, pct_counts_mt")
    adata.var["mt"] = adata.var["gene_symbols"].str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    study_metadata = study_metadata_path.read_text()
    paper_text = paper_text_path.read_text()
    quantile_table = compute_quantiles(adata)

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
    sc.pp.filter_cells(adata, min_genes=min_genes)
    if max_counts is not None:
        sc.pp.filter_cells(adata, max_counts=max_counts)
    adata = adata[adata.obs.pct_counts_mt < max_pct_mt].copy()

    stages.append(_stage_snapshot("After fixed thresholds", adata))

    print("Applying adaptive quantile trimming")
    adata.obs = adata.obs.join(quantile_table, on="latch_sample_id")

    adaptive_mask = (
        (adata.obs["n_genes_by_counts"]
             >= adata.obs[qcol("n_genes_by_counts", 0.05)]) &
        (adata.obs["n_genes_by_counts"]
             <= adata.obs[qcol("n_genes_by_counts", 0.95)]) &

        (adata.obs["total_counts"]
             >= adata.obs[qcol("total_counts", 0.05)]) &
        (adata.obs["total_counts"]
             <= adata.obs[qcol("total_counts", 0.95)]) &

        (adata.obs["pct_counts_mt"]
             <= adata.obs[qcol("pct_counts_mt", 0.95)])
    )
    print(f"Adaptive filter retains {adaptive_mask.sum()} / {adaptive_mask.size} cells")
    adata = adata[adaptive_mask, :]
    stages.append(_stage_snapshot("After adaptive thresholds", adata))
    write_anndata(adata, workdir, lcc.qc_adata_name)

    html = build_qc_report_html(stages)
    write_html_report(html, workdir, lcc.qc_report_name)

    print(f"QC pipeline finished successfully – final cell count: {adata.n_obs}")

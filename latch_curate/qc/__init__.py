from __future__ import annotations

from typing import TypedDict, Any
from pathlib import Path
import json
from textwrap import dedent
import logging

import pandas as pd
import scanpy as sc
from anndata import AnnData

from latch_curate.utils import _fig_to_base64, print_full_df_string, prompt_model

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    )

class CountsThresholds(TypedDict):
    min_genes: int
    max_total_counts: int | None
    max_pct_mt: int

def qcol(metric: str, q: float) -> str:
    return f"{metric}_q{int(q*1000):03d}"

TECH_THRESHOLDS: dict[str, CountsThresholds] = {
    # droplet‑based
    "10x Chromium v2/v3 (3′/5′)":     {"min_genes": 200,  "max_total_counts": 25_000, "max_pct_mt": 15},
    "10x Chromium Flex (3′)":         {"min_genes": 100,  "max_total_counts": 15_000, "max_pct_mt": 20},
    "Drop‑seq / Seq‑Well":            {"min_genes": 200,  "max_total_counts": 20_000, "max_pct_mt": 20},
    "inDrops":                        {"min_genes": 200,  "max_total_counts": 15_000, "max_pct_mt": 20},
    "SPLiT‑seq / sci‑RNA‑seq":        {"min_genes": 100,  "max_total_counts": 10_000, "max_pct_mt": 10},

    # nucleus droplet
    "10x Single‑nucleus (snRNA‑seq)": {"min_genes": 100,  "max_total_counts": 10_000, "max_pct_mt":  5},

    # full‑length plate
    "Smart‑seq2":                     {"min_genes": 1_000,"max_total_counts": None,   "max_pct_mt": 30},
    "Smart‑seq3":                     {"min_genes": 800,  "max_total_counts": None,   "max_pct_mt": 30},
    "FLASH‑seq":                      {"min_genes": 800,  "max_total_counts": None,   "max_pct_mt": 30},
    "CEL‑Seq2":                       {"min_genes": 500,  "max_total_counts": 50_000, "max_pct_mt": 25},

    # spatial
    "10x Visium":                     {"min_genes": 200,  "max_total_counts": 50_000, "max_pct_mt": 20},
    "Slide‑seq v2":                   {"min_genes": 100,  "max_total_counts": 12_000, "max_pct_mt": 15},
}

def _violin_plot(adata: AnnData, groupby: str | None = None):
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


def _df_to_html(df: pd.DataFrame) -> str:
    return (
        df.to_html(classes="quant", border=0, na_rep="", float_format="{:.2f}".format)
        .replace("<th>", "<th style='text-align:center'>")
        .replace("<td>", "<td style='text-align:center'>")
    )

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

def build_fixed_threshold_prompt(
    metrics: set[str],
    quantile_table: pd.DataFrame,
    study_metadata: str,
    paper_text: str,
    thresholds_by_technology: dict[str, CountsThresholds],
) -> str:

    example_output = {k: "" for k in list(metrics) + ["notes"]}

    return dedent(
        f"""
        <study_metadata>\n{study_metadata}\n</study_metadata>

        <quantile_table>\n{print_full_df_string(quantile_table)}\n</quantile_table>

        <thresholds_by_technology>\n{thresholds_by_technology}\n</thresholds_by_technology>

        I need you to construct first‑pass thresholds to QC single‑cell data.
        You get:
          • <quantile_table> per sample
          • <thresholds_by_technology> canonical suggestions
          • <study_metadata> info

        Produce ONE conservative threshold set for {metrics}.
        Output JSON exactly like {example_output}. No extra text.
        """
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


def qc_and_filter(
    adata: AnnData,
    study_metadata_path: Path,
    paper_text_path: Path,
    output_html_path: Path = Path("qc_report.html"),
) -> AnnData:
    """Runs 3‑stage QC:
       1. Raw
       2. After fixed hard thresholds
       3. After adaptive quantile trimming
    """

    logger.debug("Calculating qc metrics: mt‑genes, n_genes_by_counts, total_counts, pct_counts_mt")
    adata.var["mt"] = adata.var["gene_symbols"].str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    study_metadata = study_metadata_path.read_text()
    paper_text = paper_text_path.read_text()

    def _stage_snapshot(name: str, _adata: AnnData) -> dict[str, Any]:
        return {
            "name": name,
            "total": _adata.n_obs,
            "counts": _adata.obs["latch_sample_id"].value_counts().to_dict(),
            "global_fig": _fig_to_base64(_violin_plot(_adata)),
            "sample_fig": _fig_to_base64(_violin_plot(_adata, groupby="latch_sample_id")),
            "qtable": _df_to_html(compute_quantiles(_adata)),
        }

    stages: list[dict[str, Any]] = []
    stages.append(_stage_snapshot("Before filtering", adata))

    quantile_table = compute_quantiles(adata)
    prompt = build_fixed_threshold_prompt(
        {"min_genes", "max_total_counts", "max_pct_mt"},
        quantile_table,
        study_metadata,
        paper_text,
        TECH_THRESHOLDS,
    )

    logger.info("Requesting fixed thresholds from language model")
    while True:
        try:
            fixed = json.loads(prompt_model(prompt))
            min_genes = int(fixed["min_genes"])
            max_counts = None if fixed["max_total_counts"] is None else int(fixed["max_total_counts"])
            max_pct_mt = int(fixed["max_pct_mt"])
            break
        except Exception:
            print("[filter] invalid model response – retrying …")

    logger.debug("Applying fixed thresholds")
    sc.pp.filter_cells(adata, min_genes=min_genes)
    if max_counts is not None:
        sc.pp.filter_cells(adata, max_counts=max_counts)
    adata = adata[adata.obs.pct_counts_mt < max_pct_mt].copy()

    stages.append(_stage_snapshot("After fixed thresholds", adata))

    logger.info("Applying adaptive quantile trimming")
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
    logger.debug("Adaptive filter retains %d / %d cells", adaptive_mask.sum(),
                 adaptive_mask.size)
    adata = adata[adaptive_mask, :]

    stages.append(_stage_snapshot("After adaptive thresholds", adata))

    html = build_qc_report_html(stages)
    output_html_path.write_text(html)

    print(f"[qc] report written to {output_html_path.resolve()}")
    logger.info("QC pipeline finished successfully – final cell count: %d", adata.n_obs)
    return output_html_path

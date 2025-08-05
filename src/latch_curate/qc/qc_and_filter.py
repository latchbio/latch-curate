from __future__ import annotations

import warnings
from typing import Any, TypedDict
from pathlib import Path
from textwrap import dedent
import json
import yaml

import pandas as pd
import scanpy as sc
from anndata import AnnData, ImplicitModificationWarning
from matplotlib import pyplot as plt

from latch_curate.utils import _fig_to_base64, write_anndata, write_html_report, df_to_str
from latch_curate.constants import latch_curate_constants as lcc
from latch_curate.utils import _df_to_html as _df_to_html_base
from latch_curate.tinyrequests import post
from latch_curate.config import user_config

def _df_to_html(df: pd.DataFrame, class_name: str = "quant") -> str:
    html = _df_to_html_base(df)
    if class_name != "quant":
        html = html.replace('class="quant"', f'class="{class_name}"')
    return html

statistics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]

class StatisticInterval(TypedDict):
    statistic_name: str
    interval:  tuple[float | None, float | None]
    reasoning: str

def qcol(metric: str, q: float) -> str:
    return f"{metric}_q{int(q*1000):03d}"

def _violin_plot(adata: AnnData, groupby: str | None = None):
    if adata.n_obs == 0 or (groupby and adata.obs[groupby].nunique() == 0):
        return plt.figure()
    if groupby and pd.api.types.is_categorical_dtype(adata.obs[groupby]):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=ImplicitModificationWarning)
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

def _violin_plot_single(adata: AnnData, metric: str, groupby: str | None = None, thresholds: dict | None = None):
    if adata.n_obs == 0 or (groupby and adata.obs[groupby].nunique() == 0):
        return plt.figure()
    if groupby and pd.api.types.is_categorical_dtype(adata.obs[groupby]):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=ImplicitModificationWarning)
            adata.obs[groupby] = adata.obs[groupby].cat.remove_unused_categories()

    fig, ax = plt.subplots(figsize=(8, 5))
    sc.pl.violin(adata, metric, groupby=groupby, ax=ax, show=False, jitter=0.4)

    if thresholds:
        for sample_id, bounds in thresholds.items():
            if bounds.get('lower') is not None:
                ax.axhline(y=bounds['lower'], color='red', linestyle='--', alpha=0.7, label=f'{sample_id} lower')
            if bounds.get('upper') is not None:
                ax.axhline(y=bounds['upper'], color='red', linestyle='--', alpha=0.7, label=f'{sample_id} upper')

    ax.set_title(f'{metric} distribution')
    plt.tight_layout()
    return fig


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

def compute_quantiles_by_metric(adata: AnnData) -> dict[str, pd.DataFrame]:
    tail_quants = {"low": [0.01, 0.05], "high": [0.90, 0.95, 0.99, 0.995]}
    metrics = {
        "n_genes_by_counts": ["low", "high"],
        "total_counts":      ["low", "high"],
        "pct_counts_mt":     ["high"],
    }

    result = {}
    for metric, tails in metrics.items():
        metric_quantiles = {}
        for t in tails:
            for q in tail_quants[t]:
                col_name = f"q{int(q*1000):03d}"
                metric_quantiles[col_name] = adata.obs.groupby("latch_sample_id")[metric].quantile(q)
        result[metric] = pd.DataFrame(metric_quantiles)

    return result

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
        stage_html = [
            f"<h2>{s['name']}</h2>",
            "<div class='summary'>",
            f"  <p><strong>Total cells:</strong> {s['total']:,}</p>",
            "  <table class='sample-counts'><tr><th>Sample ID</th><th>Cell Count</th></tr>",
        ]
        stage_html.extend(f"    <tr><td>{sid}</td><td>{c:,}</td></tr>" for sid, c in s["counts"].items())
        stage_html.append("  </table></div>")

        if s.get('thresholds'):
            stage_html.append("<div class='thresholds'>")
            stage_html.append("<h3>Applied Thresholds</h3>")
            stage_html.append(s['thresholds'])
            stage_html.append("</div>")

        if s.get('metric_sections'):
            stage_html.append("<div class='metrics'>")
            stage_html.extend(s['metric_sections'])
            stage_html.append("</div>")
        else:
            stage_html.extend([
                "<h3>Quantiles</h3>", s.get("qtable", ""),
                "<h3>Distributions</h3>",
                f"<img src='data:image/png;base64,{s['global_fig']}' alt='{s['name']} global'/>",
                f"<img src='data:image/png;base64,{s['sample_fig']}' alt='{s['name']} by-sample'/>",
            ])

        if s.get("interval_html"):
            stage_html.extend([
                "<h3>Adaptive intervals &amp; reasoning</h3>",
                s["interval_html"]
            ])

        rows_html.append("\n".join(stage_html))

    body = "\n".join(rows_html)
    return dedent(
        f"""<!DOCTYPE html>
        <html><head><meta charset='utf-8'><title>{title}</title>
        <style>
          body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; 
                  max-width: 1200px; margin: 0 auto; padding: 20px; 
                  background-color: #f8f9fa; }}
          h1 {{ color: #2c3e50; text-align: center; margin-bottom: 40px; }}
          h2 {{ color: #34495e; border-bottom: 2px solid #3498db; padding-bottom: 10px; 
                margin-top: 40px; }}
          h3 {{ color: #34495e; margin-top: 30px; }}
          .summary {{ background: white; padding: 20px; border-radius: 8px; 
                      box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 20px; }}
          .metrics {{ display: flex; flex-direction: column; gap: 30px; }}
          .metric-section {{ background: white; padding: 20px; border-radius: 8px; 
                             box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
          .metric-header {{ display: flex; align-items: center; gap: 20px; 
                           border-bottom: 1px solid #e0e0e0; padding-bottom: 15px; 
                           margin-bottom: 20px; }}
          .metric-name {{ font-size: 1.3em; font-weight: bold; color: #2c3e50; }}
          .metric-content {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
          table {{ border-collapse: collapse; width: 100%; background: white; }}
          th, td {{ border: 1px solid #e0e0e0; padding: 8px 12px; text-align: center; }}
          th {{ background-color: #ecf0f1; font-weight: bold; color: #2c3e50; }}
          tr:nth-child(even) {{ background-color: #f8f9fa; }}
          .sample-counts {{ max-width: 400px; margin: 20px auto; }}
          .quantile-table {{ font-size: 0.9em; }}
          .threshold-line {{ color: #e74c3c; font-weight: bold; }}
          img {{ max-width: 100%; height: auto; border-radius: 4px; 
                 box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
          .thresholds {{ background: #fff3cd; padding: 15px; border-radius: 8px; 
                         border: 1px solid #ffeeba; margin-bottom: 20px; }}
          .thresholds h3 {{ margin-top: 0; color: #856404; }}
          .threshold-item {{ margin: 5px 0; }}
          .adaptive-thresholds {{ margin-top: 20px; }}
          .threshold-stat {{ background: white; padding: 20px; border-radius: 8px; 
                             box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 20px; }}
          .threshold-stat h4 {{ margin-top: 0; color: #2c3e50; border-bottom: 1px solid #e0e0e0; 
                                padding-bottom: 10px; }}
          .threshold-table {{ font-size: 0.95em; }}
          .bound-value {{ font-family: monospace; font-weight: bold; color: #e74c3c; }}
          .reasoning {{ font-style: italic; color: #7f8c8d; }}
          .quantile-panel, .plot-panel {{ background: #f8f9fa; padding: 15px; border-radius: 6px; }}
          .quantile-panel h4, .plot-panel h4 {{ margin-top: 0; color: #34495e; }}
          @media (max-width: 768px) {{
            .metric-content {{ grid-template-columns: 1fr; }}
          }}
        </style></head><body>
          <h1>{title}</h1>
          {body}
        </body></html>"""
    )

def _intervals_to_html(intervals: dict[str, list[StatisticInterval]]) -> str:
    if not intervals:
        return ""

    by_stat = {}
    for sample, stats in intervals.items():
        for stat in stats:
            if stat["statistic_name"] not in by_stat:
                by_stat[stat["statistic_name"]] = {}
            by_stat[stat["statistic_name"]][sample] = stat

    metric_names = {
        "n_genes_by_counts": "Number of Genes per Cell",
        "total_counts": "Total UMI Counts per Cell", 
        "pct_counts_mt": "Mitochondrial Gene Percentage"
    }

    rows = ['<div class="adaptive-thresholds">']

    for stat_name, samples in by_stat.items():
        display_name = metric_names.get(stat_name, stat_name)
        rows.append('<div class="threshold-stat">')
        rows.append(f'<h4>{display_name}</h4>')
        rows.append('<table class="threshold-table">')
        rows.append('<tr><th>Sample</th><th>Lower Bound</th><th>Upper Bound</th><th>Reasoning</th></tr>')

        for sample, stat in samples.items():
            lower = stat["interval"][0]
            upper = stat["interval"][1]
            lower_str = f'{lower:.2f}' if lower is not None else '-'
            upper_str = f'{upper:.2f}' if upper is not None else '-'

            rows.append(
                f'<tr>'
                f'<td><strong>{sample}</strong></td>'
                f'<td class="bound-value">{lower_str}</td>'
                f'<td class="bound-value">{upper_str}</td>'
                f'<td class="reasoning">{stat["reasoning"]}</td>'
                f'</tr>'
            )

        rows.append('</table>')
        rows.append('</div>')

    rows.append('</div>')
    return '\n'.join(rows)

def _stage_snapshot(
    name: str,
    _adata: AnnData,
    *,
    intervals: dict[str, list[StatisticInterval]] | None = None,
    fixed_thresholds: dict | None = None,
) -> dict[str, Any]:
    threshold_data = {}
    if intervals:
        for sample_id, stats in intervals.items():
            for stat in stats:
                if stat["statistic_name"] not in threshold_data:
                    threshold_data[stat["statistic_name"]] = {}
                threshold_data[stat["statistic_name"]][sample_id] = {
                    'lower': stat["interval"][0],
                    'upper': stat["interval"][1]
                }

    metric_sections = []
    quantiles_by_metric = compute_quantiles_by_metric(_adata)

    metric_names = {
        "n_genes_by_counts": "Number of Genes per Cell",
        "total_counts": "Total UMI Counts per Cell",
        "pct_counts_mt": "Mitochondrial Gene Percentage"
    }

    for metric, display_name in metric_names.items():
        section_html = [
            '<div class="metric-section">',
            '  <div class="metric-header">',
            f'    <div class="metric-name">{display_name}</div>',
            '  </div>',
            '  <div class="metric-content">',
            '    <div class="quantile-panel">',
            '      <h4>Quantiles by Sample</h4>',
            f'      {_df_to_html(quantiles_by_metric[metric], class_name="quantile-table")}',
            '    </div>',
            '    <div class="plot-panel">',
            '      <h4>Distribution</h4>',
        ]

        plot_thresholds = threshold_data.get(metric)
        fig = _violin_plot_single(_adata, metric, groupby="latch_sample_id", thresholds=plot_thresholds)
        section_html.append(f'      <img src="data:image/png;base64,{_fig_to_base64(fig)}" alt="{metric} distribution"/>')

        section_html.extend([
            '    </div>',
            '  </div>',
            '</div>'
        ])

        metric_sections.append('\n'.join(section_html))

    threshold_html = ""
    if fixed_thresholds:
        threshold_items = []
        if fixed_thresholds.get('min_genes'):
            threshold_items.append(f'<div class="threshold-item">Minimum genes: {fixed_thresholds["min_genes"]}</div>')
        if fixed_thresholds.get('max_counts'):
            threshold_items.append(f'<div class="threshold-item">Maximum counts: {fixed_thresholds["max_counts"]:,}</div>')
        if fixed_thresholds.get('max_pct_mt'):
            threshold_items.append(f'<div class="threshold-item">Maximum MT%: {fixed_thresholds["max_pct_mt"]}%</div>')
        threshold_html = '\n'.join(threshold_items)

    return {
        "name": name,
        "total": _adata.n_obs,
        "counts": _adata.obs["latch_sample_id"].value_counts().to_dict(),
        "metric_sections": metric_sections,
        "thresholds": threshold_html,
        "interval_html": _intervals_to_html(intervals or {}),
        "global_fig": _fig_to_base64(_violin_plot(_adata)),
        "sample_fig": _fig_to_base64(_violin_plot(_adata, groupby="latch_sample_id")),
        "qtable": _df_to_html(compute_quantiles(_adata)),
    }


def qc_and_filter(
    adata: AnnData,
    study_metadata_path: Path,
    paper_text_path: Path,
    workdir: Path,
    use_params: bool
):
    params_path = workdir / lcc.qc_params_name
    if use_params:
        assert params_path.exists(), f"Missing parameters file: {params_path}"
        with open(params_path) as f:
            params = yaml.safe_load(f)

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
    post_fixed_quantile_table = compute_quantiles(adata)

    fixed_threshold_dict = {
        'min_genes': min_genes,
        'max_counts': max_counts,
        'max_pct_mt': max_pct_mt
    }

    stages.append(_stage_snapshot("After fixed thresholds", adata, fixed_thresholds=fixed_threshold_dict))

    print("Applying adaptive thresholds")

    interval_data: dict[str, list[StatisticInterval]] = {}

    if use_params:
        interval_data = params['adaptive']
        print("Using existing adaptive thresholds from params file")
    else:
        print("Requesting adaptive thresholds from model")
        resp = post(
            f"{lcc.nucleus_url}/{lcc.get_adaptive_qc_thresholds_endpoint}",
            {
                "statistics": statistics,
                "quantile_table": df_to_str(post_fixed_quantile_table),
                "sample_list": list(set(adata.obs['latch_sample_id'])),
                "session_id": -1
            },
            headers = {"Authorization": f"Latch-SDK-Token {user_config.token}"}
        )
        try:
            data = resp.json()['data']['interval_data']
            for sample_name, x in data.items():
                interval_data[sample_name] = []
                for v in x:
                    assert "statistic_name" in v
                    assert "interval" in v
                    assert "reasoning" in v
                    stat_data = StatisticInterval(statistic_name=v["statistic_name"], interval=v["interval"], reasoning=v["reasoning"])
                    interval_data[sample_name].append(stat_data)
        except Exception as e:
            print(e.with_traceback())
            print(f"Invalid model response: {e}. Trying again")

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
        with open(params_path, "w") as f:
            yaml.safe_dump(param_data, f, default_flow_style=False, indent=2)
        print(f"QC parameters written to {params_path}")
    print(f"QC pipeline finished successfully – final cell count: {adata.n_obs}")

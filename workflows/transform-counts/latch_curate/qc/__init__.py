from typing import Dict

from pathlib import Path
from anndata import AnnData
import scanpy as sc
import matplotlib.pyplot as plt
from textwrap import dedent

from utils import _fig_to_base64

max_pct_counts_mt = 15
max_counts_per_cell = 20_000
min_genes = 100

stringent_min_genes = 250

# note(kenny): There is a lot to configure and will parameterize over time based on need.
# Future metrics:
# - % ribosomal
# - ambient signatures
# - doublet score
# - cell-cycle phase
# - complexity ratio : (n_genes / total_counts)
# - novelity : 1 - log(n_genes) log(total_counts)
# Potential directions: filtering by adaptive thresholds computed from each

# use adaptive thresholds
# qs = adata.obs.groupby('latch_sample_id').agg({
#     'n_genes_by_counts':   [lambda x: x.quantile(0.05)],
#     'total_counts':        [lambda x: x.quantile(0.95)],
#     'pct_counts_mt':       [lambda x: x.quantile(0.95)],
# })
# 
# qs.columns = ['min_genes_q','max_counts_q','max_mt_q']
# qs = qs.reset_index()

# adata.obs = adata.obs.merge(qs, on='latch_sample_id', how='left')
# 
# mask = (
#     (adata.obs['n_genes_by_counts'] >= adata.obs['min_genes_q']) &
#     (adata.obs['total_counts']      <= adata.obs['max_counts_q']) &
#     (adata.obs['pct_counts_mt']     <= adata.obs['max_mt_q'])
# )

# adata_q = adata[mask].copy()


def build_qc_report_html(
    total_before: int,
    counts_before: Dict[str, int],
    total_after: int,
    counts_after: Dict[str, int],
    figs: Dict[str, plt.Figure],
    title: str = "QC & Filtering Report"
) -> str:
    imgs = {name: _fig_to_base64(fig) for name, fig in figs.items()}

    def _rows(counts: Dict[str,int]) -> str:
        return "".join(f"<tr><td>{s}</td><td>{c}</td></tr>"
                       for s, c in counts.items())

    return dedent(
        f"""<!DOCTYPE html>
            <html>
            <head>
              <meta charset="utf-8">
              <title>{title}</title>
              <style>
                body {{ font-family: sans-serif; max-width: 800px; margin: auto; }}
                h2 {{ border-bottom: 1px solid #ccc; }}
                table {{ border-collapse: collapse; width: 100%; margin-bottom: 1em; }}
                th, td {{ border: 1px solid #ddd; padding: 4px 8px; text-align: center; }}
                img {{ max-width: 100%; height: auto; margin-bottom: 1em; }}
              </style>
            </head>
            <body>
              <h1>{title}</h1>

              <h2>Before Filtering</h2>
              <p><strong>Total cells:</strong> {total_before}</p>
              <h3>Per-sample counts</h3>
              <table>
                <tr><th>Sample ID</th><th>Count</th></tr>
                {_rows(counts_before)}
              </table>
              <h3>Distributions</h3>
              <img src="data:image/png;base64,{imgs['before_global']}" alt="Before: global violins"/>
              <img src="data:image/png;base64,{imgs['before_by_sample']}" alt="Before: by-sample violins"/>

              <h2>After Filtering</h2>
              <p><strong>Total cells:</strong> {total_after}</p>
              <h3>Per-sample counts</h3>
              <table>
                <tr><th>Sample ID</th><th>Count</th></tr>
                {_rows(counts_after)}
              </table>
              <h3>Distributions</h3>
              <img src="data:image/png;base64,{imgs['after_global']}" alt="After: global violins"/>
              <img src="data:image/png;base64,{imgs['after_by_sample']}" alt="After: by-sample violins"/>
            </body>
            </html>
        """)

def construct_violin_plot(adata: AnnData):
    grid = sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, show=False)
    print(grid)
    return grid.fig

def construct_per_sample_violin_plot(adata: AnnData):
    grid = sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, groupby='latch_sample_id', rotation=45, show=False)
    return grid[0].figure

def qc_and_filter(adata: AnnData, output_html: Path = Path("qc_report.html")):

    print("[qc] calculating metrics")
    adata.var["mt"] = adata.var["gene_symbols"].str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    return

    total_before = adata.n_obs
    vc_before = adata.obs['latch_sample_id'].value_counts()

    fig_bg = construct_violin_plot(adata)
    fig_bs = construct_per_sample_violin_plot(adata)

    print("[filter] first pass filter with fixed thresholds")
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, max_counts=max_counts_per_cell)
    adata = adata[adata.obs.pct_counts_mt < max_pct_counts_mt, :] 

    counts_by_sample = adata.obs['latch_sample_id'].value_counts().to_dict()
    samples_to_filter = [s for s,c in counts_by_sample.items() if c >= 7_000]

    if len(samples_to_filter) > 0:
        print("[filter] stringent filter on outlier samples")
        mask = (~adata.obs['latch_sample_id'].isin(samples_to_filter)) | (adata.obs['n_genes_by_counts'] > stringent_min_genes)
        adata = adata[mask, :]

    fig_ag = construct_violin_plot(adata)
    fig_as = construct_per_sample_violin_plot(adata)

    total_after = adata.n_obs
    vc_after = adata.obs['latch_sample_id'].value_counts()

    html = build_qc_report_html(
        total_before=total_before,
        counts_before=vc_before.to_dict(),
        total_after=total_after,
        counts_after=vc_after.to_dict(),
        figs={
            'before_global': fig_bg,
            'before_by_sample': fig_bs,
            'after_global': fig_ag,
            'after_by_sample': fig_as
        }
    )
    Path(output_html).write_text(html)

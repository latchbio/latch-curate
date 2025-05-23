import scanpy as sc
from anndata import AnnData
import scanpy.external as sce
import matplotlib.pyplot as plt

from pathlib import Path
from textwrap import dedent

from latch_curate.utils import _fig_to_base64


def build_transform_report_html(
    figs: dict[str, plt.Figure],
    output_path: Path,
    title: str = "Transform & Integration Report"
) -> None:
    """Embed a set of named figures into a single HTML page."""
    imgs = {name: _fig_to_base64(fig) for name, fig in figs.items()}

    body = "\n".join(
        f"<h2>{name}</h2>\n"
        f'<img src="data:image/png;base64,{data}" alt="{name}"/>'
        for name, data in imgs.items()
    )

    html = dedent(f"""\
        <!DOCTYPE html>
        <html>
        <head>
          <meta charset="utf-8">
          <title>{title}</title>
          <style>
            body {{ font-family: sans-serif; max-width: 800px; margin: auto; }}
            h2 {{ margin-top: 1.5em; border-bottom: 1px solid #ccc; }}
            img {{ width: 100%; height: auto; margin-bottom: 1em; }}
          </style>
        </head>
        <body>
          <h1>{title}</h1>
          {body}
        </body>
        </html>
    """)

    output_path.write_text(html)

def transform_counts(adata: AnnData, use_scrublet: bool = True, output_html_path: Path = Path("transform_report.html")):
    if use_scrublet:
        print("[scrublet] starting scrublet")
        sc.pp.scrublet(adata)
        print("[scrublet] finished scrublet")

    print("[pca] starting counts transform")
    adata.raw = adata.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    print("[pca] finished counts transform")

    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="latch_sample_id")
    hvg_fig = sc.pl.highly_variable_genes(adata, show=False)

    print("[pca] starting pca")
    sc.tl.pca(adata)
    # sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)

    pca_fig  = sc.pl.pca(
        adata,
        color=["latch_sample_id", "latch_sample_id", "total_counts", "total_counts"],
        dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
        ncols=2,
        size=2,
        show=False
    )
    print("[pca] finished pca")

    print("[integrate] starting harmony")
    sce.pp.harmony_integrate(adata, 'latch_sample_id')
    harmony_pca_fig = sc.pl.embedding(
        adata, 
        basis='X_pca_harmony', 
        color=['latch_sample_id', 'total_counts'],
        title=['Harmony PCA: Sample', 'Harmony PCA: Total Counts'],
        show=False
    )
    pca_fig = sc.pl.embedding(
        adata, 
        basis='X_pca', 
        color=['latch_sample_id', 'total_counts'],
        title=['PCA: Sample', 'PCA: Total Counts'],
        show=False
    )
    adata.obsm["X_pca"] = adata.obsm["X_pca_harmony"]
    del adata.obsm["X_pca_harmony"]
    print("[integrate] finished harmony")

    print("[cluster] starting nearest neighbors")
    sc.pp.neighbors(adata)
    print("[cluster] finished nearest neighbors")

    print("[dim reduction] starting umap")
    sc.tl.umap(adata)
    umap_fig = sc.pl.umap(
        adata,
        color="latch_sample_id",
        size=2,
        show=False
    )
    print("[dim reduction] finished umap")

    leiden_resolutions = [0.50, 1.00, 2.00]
    leiden_figs = {}
    for res in leiden_resolutions:
       key = f"leiden_res_{res:4.2f}"
       sc.tl.leiden(adata, key_added=key, resolution=res, flavor="igraph")
       leiden_figs[f"Leiden {res:.2f}"]= sc.pl.umap(
           adata,
           color=[key],
           legend_loc="on data",
           show=False
       )

    build_transform_report_html(
        {
               "Highly Variable Genes": hvg_fig,
               "PCA":                   pca_fig,
               "Harmony PCA":           harmony_pca_fig,
               "UMAP":                  umap_fig,
               **leiden_figs,
        },
        output_html_path
    )

    return adata, output_html_path

from pathlib import Path
from textwrap import dedent

import matplotlib.pyplot as plt
import scanpy as sc
import scanpy.external as sce
from anndata import AnnData

from latch_curate.utils import _fig_to_base64, write_html_report, write_anndata
from latch_curate.cell_typing.marker_genes import marker_genes, remove_absent_symbols
from latch_curate.constants import latch_curate_constants as lcc

def build_transform_report_html(
    figs: dict[str, plt.Figure],
    title: str = "Transform & Integration Report",
) -> str:
    imgs = {name: _fig_to_base64(fig) for name, fig in figs.items()}
    body = "\n".join(
        f"<h2>{name}</h2>\n"
        f'<img src="data:image/png;base64,{data}" alt="{name}"/>'
        for name, data in imgs.items()
    )
    html = dedent(
        f"""<!DOCTYPE html>
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
        </html>"""
    )
    return html


def transform_counts(
    adata: AnnData,
    workdir: Path,
    use_scrublet: bool = True,
):
    workdir.mkdir(exist_ok=True)

    if use_scrublet:
        print("starting scrublet")
        sc.pp.scrublet(adata)
        print("finished scrublet")

    print("normalizing counts")
    adata.raw = adata.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    print("finished normalization")

    print("finding highly variable genes")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="latch_sample_id")
    hvg_fig = sc.pl.highly_variable_genes(adata, show=False)

    print("running PCA")
    sc.tl.pca(adata)
    pca_fig = sc.pl.embedding(
        adata,
        basis="X_pca",
        color=["latch_sample_id", "total_counts"],
        title=["PCA: Sample", "PCA: Total Counts"],
        show=False,
    )
    print("finished PCA")

    print("running Harmony integration")
    sce.pp.harmony_integrate(adata, "latch_sample_id")
    harmony_pca_fig = sc.pl.embedding(
        adata,
        basis="X_pca_harmony",
        color=["latch_sample_id", "total_counts"],
        title=["Harmony PCA: Sample", "Harmony PCA: Total Counts"],
        show=False,
    )
    adata.obsm["X_pca"] = adata.obsm.pop("X_pca_harmony")
    print("finished Harmony integration")

    print("computing neighbors")
    sc.pp.neighbors(adata)
    print("finished neighbors")

    print("computing UMAP")
    sc.tl.umap(adata)
    umap_fig = sc.pl.umap(adata, color="latch_sample_id", size=2, show=False)
    print("finished UMAP")

    print("generating Leiden clusterings")
    leiden_resolutions = [0.50, 1.00, 2.00]
    leiden_figs: dict[str, plt.Figure] = {}
    for res in leiden_resolutions:
        key = f"leiden_res_{res:4.2f}"
        sc.tl.leiden(adata, key_added=key, resolution=res, flavor="igraph")
        leiden_figs[f"Leiden {res:.2f}"] = sc.pl.umap(
            adata,
            color=[key],
            legend_loc="on data",
            show=False,
        )
    print("finished Leiden clusterings")

    print("generating marker-gene dotplot")

    filtered_markers = remove_absent_symbols(adata, marker_genes)
    dp = sc.pl.dotplot(
        adata,
        filtered_markers,
        gene_symbols="gene_symbols",
        groupby="leiden_res_0.50",
        standard_scale="var",
        show=False,
    )
    dotplot_fig = next(iter(dp.values())).figure
    print("finished dotplot")


    html = build_transform_report_html(
        {
            "Highly Variable Genes": hvg_fig,
            "PCA": pca_fig,
            "Harmony PCA": harmony_pca_fig,
            "UMAP": umap_fig,
            "Marker Gene DotPlot": dotplot_fig,
            **leiden_figs,
        },
    )
    write_html_report(html, workdir, lcc.transform_report_name)
    write_anndata(adata, workdir, lcc.transform_adata_name)

    print("transform pipeline complete")

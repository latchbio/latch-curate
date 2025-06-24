import json
from pathlib import Path
from textwrap import dedent

import matplotlib.pyplot as plt
import scanpy as sc
from anndata import AnnData

from latch_curate.cell_typing.marker_genes import marker_genes, remove_absent_symbols
from latch_curate.cell_typing.vocab import cell_typing_vocab
from latch_curate.config import user_config
from latch_curate.utils import _fig_to_base64, write_html_report, write_anndata
from latch_curate.tinyrequests import post
from latch_curate.constants import latch_curate_constants as lcc

cluster_key = "leiden_res_0.50"

def construct_diff_exp_json(adata: AnnData, top_n: int = 20) -> dict:
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method="wilcoxon", n_genes=50)
    res = adata.uns["rank_genes_groups"]
    clusters = res["names"].dtype.names
    out: dict[str, list[dict]] = {}
    for c in clusters:
        genes = res["names"][c][:top_n]
        lfc = res["logfoldchanges"][c][:top_n]
        p = res["pvals_adj"][c][:top_n]
        out[c] = [
            {
                "gene": adata.var.loc[g, "gene_symbols"],
                "logfoldchange": float(lfc[i]),
                "pval_adj": float(p[i]),
            }
            for i, g in enumerate(genes)
        ]
    return out


def _build_report_html(ann: dict, reasoning_md: str, dotplot_fig: plt.Figure,
                       umap_cell_type_fig: plt.Figure, umap_leiden_fig: plt.Figure) -> None:
    dotplot_img = _fig_to_base64(dotplot_fig)
    umap_cell_type_img = _fig_to_base64(umap_cell_type_fig)
    umap_leiden_img = _fig_to_base64(umap_leiden_fig)
    html = dedent(
        f"""<!DOCTYPE html>
        <html><head><meta charset="utf-8"><title>Cell-typing Report</title>
        <style>
            body {{font-family:sans-serif;max-width:900px;margin:auto}}
            pre {{white-space:pre-wrap;border:1px solid #ddd;padding:1em;background:#fafafa}}
            img {{max-width:100%;height:auto;margin:1em 0}}
        </style></head><body>
          <h1>Cell-typing Report</h1>

          <h2>Annotations</h2>
          <pre>{json.dumps(ann, indent=2)}</pre>

          <h2>Reasoning</h2>
          <pre>{reasoning_md}</pre>

          <h2>Marker Dot Plot</h2>
          <img src="data:image/png;base64,{dotplot_img}" alt="DotPlot"/>

          <h2>UMAPs</h2>
          <img src="data:image/png;base64,{umap_cell_type_img}" alt="DotPlot"/>
          <img src="data:image/png;base64,{umap_leiden_img}" alt="DotPlot"/>
        </body></html>"""
    )
    return html


def type_cells(
    adata: AnnData,
    workdir: Path,
):
    workdir.mkdir(exist_ok=True)

    print("starting differential expression")
    de_json = construct_diff_exp_json(adata)

    print("Requesting cell types from model")
    resp = post(
        f"{lcc.nucleus_url}/{lcc.get_cell_types_endpoint}",
        {
            "marker_genes": json.dumps(marker_genes),
            "de_json": json.dumps(de_json),
            "cell_typing_vocab": json.dumps(cell_typing_vocab),
            "metadata": json.dumps({"step": "type-cells", "project":
                                    workdir.resolve().parent.name}),
            "session_id": -1
        },
        headers = {"Authorization": f"Latch-SDK-Token {user_config.token}"}
    )
    try:
        data = resp.json()['data']
        annotations = data['annotations']
        reasoning = data['reasoning']
    except KeyError:
        raise ValueError(f'Malformed response data: {resp.json()}')

    print("mapping annotations onto AnnData")
    adata.obs["latch_cell_type_lvl_1"] = adata.obs[cluster_key].astype(str).map(annotations)

    print("creating marker dotplot")
    filtered_markers = remove_absent_symbols(adata, marker_genes)
    dp = sc.pl.dotplot(
        adata,
        filtered_markers,
        gene_symbols="gene_symbols",
        groupby=cluster_key,
        standard_scale="var",
        show=False,
    )
    dotplot_fig = next(iter(dp.values())).figure
    print("finished dotplot")

    print("creating UMAP coloured by cell-type annotations")
    umap_cell_type_fig = sc.pl.umap(
        adata,
        color="latch_cell_type_lvl_1",
        size=2,
        legend_loc="on data",
        show=False,
    )
    print("finished UMAP")

    print("creating UMAP coloured by leiden annotations")
    umap_leiden_fig = sc.pl.umap(
        adata,
        color=cluster_key,
        size=2,
        legend_loc="on data",
        show=False,
    )
    print("finished UMAP")

    html = _build_report_html(annotations, reasoning, dotplot_fig, umap_cell_type_fig, umap_leiden_fig)

    write_html_report(html, workdir, lcc.type_cells_report_name)
    write_anndata(adata, workdir, lcc.type_cells_adata_name)
    print("cell-typing pipeline complete")

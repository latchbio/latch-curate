import logging
import json
from pathlib import Path
from textwrap import dedent

import matplotlib.pyplot as plt
import scanpy as sc
from anndata import AnnData

from latch_curate.cell_typing.marker_genes import marker_genes, remove_absent_symbols
from latch_curate.cell_typing.vocab import cell_typing_vocab
from latch_curate.utils import _fig_to_base64, prompt_model

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

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


def build_de_prompt(diff_exp: dict, markers: dict, vocab: dict) -> str:
    example = {
        "annotations": {"0": "endothelial cell/CL:0000115"},
        "reasoning": "### Cluster 0\n…",
    }
    return dedent(
        f"""
        <markers>
        {markers}
        </markers>

        <vocab>
        {vocab}
        </vocab>

        <diff-exp-data>
        {diff_exp}
        </diff-exp-data>

        Annotate each cluster in <diff-exp-data> with a cell type chosen **exactly** from <vocab>.
        • Treat clusters that only show S100A8/9 epithelial expression as epithelial.
        • Treat clusters with NK-cell markers as T cells.
        • If uncertain, use `unknown/`.

        Return JSON with keys `annotations` and `reasoning`.
        `reasoning` must be a **full markdown document with newlines**, like the example below.

        <example>
        {json.dumps(example, indent=2)}
        </example>
        """
    )


def _build_report_html(ann: dict, reasoning_md: str, dotplot_fig: plt.Figure, umap_fig: plt.Figure, out: Path) -> None:
    dotplot_img = _fig_to_base64(dotplot_fig)
    umap_img = _fig_to_base64(umap_fig)
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

          <h2>Marker Gene UMAP</h2>
          <img src="data:image/png;base64,{umap_img}" alt="DotPlot"/>
        </body></html>"""
    )
    out.write_text(html)
    logger.info("wrote cell-typing report to %s", out.resolve())


def type_cells(
    adata: AnnData,
    output_html_path: Path = Path("cell_typing_report.html"),
) -> Path:
    logger.info("starting differential expression")
    de_json = construct_diff_exp_json(adata)
    prompt = build_de_prompt(de_json, marker_genes, cell_typing_vocab)

    logger.info("sending prompt to LLM")
    while True:
        try:
            resp = json.loads(prompt_model(prompt))
            annotations: dict = resp["annotations"]
            reasoning: str = resp["reasoning"]
            break
        except Exception:
            logger.warning("invalid model response, retrying…")

    logger.info("mapping annotations onto AnnData")
    adata.obs["latch_cell_type_lvl_1"] = adata.obs[cluster_key].astype(str).map(annotations)

    logger.info("creating marker dotplot")
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
    logger.info("finished dotplot")

    logger.info("creating UMAP coloured by cell-type annotations")
    umap_fig = sc.pl.umap(
        adata,
        color="latch_cell_type_lvl_1",
        size=2,
        legend_loc="on data",
        show=False,
    )
    logger.info("finished UMAP")

    _build_report_html(annotations, reasoning, dotplot_fig, umap_fig, output_html_path)
    logger.info("cell-typing pipeline complete")
    return output_html_path

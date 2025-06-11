from pathlib import Path
from typing import Literal

import scanpy as sc
import click
import pandas as pd

from latch.ldata._transfer.node import get_node_data

from latch_curate.constants import latch_curate_constants as lcc
from latch_curate.download import construct_study_metadata, download_gse_supps, get_subseries_ids
from latch_curate.construct import construct_counts as _construct_counts
from latch_curate.qc import qc_and_filter
from latch_curate.transform import transform_counts
from latch_curate.cell_typing import type_cells as _type_cells
from latch_curate.harmonize import harmonize_metadata as _harmonize_metadata
from latch_curate.lint import lint_anndata
from latch_curate.publish.build import build_publish_data
from latch_curate.publish.queries import build_publish_queries

StepwiseAction = Literal["run", "validate"]
stepwise_actions: list[StepwiseAction] = ["run", "validate"]

@click.group("latch-curate", context_settings={"max_content_width": 160},
             invoke_without_command=True)
@click.version_option(package_name=lcc.pkg_name)
def main():
    """Tools to curate single cell sequencing data. """
    # local_ver = parse(get_local_package_version())
    # latest_ver = parse(get_latest_package_version())
    # if local_ver < latest_ver:
    #     click.secho(
    #         dedent(f"""
    #             WARN: Your local version of latch-curate ({local_ver}) is out of date. This may result in unexpected behavior.
    #             Please upgrade to the latest version ({latest_ver}) using `python3 -m pip install --upgrade latch-curate`.
    #             """).strip("\n"),
    #         fg="yellow",
    #     )

    # crash_handler.init()

# @click.command("init")
# def init(gse_id: str):
#     ...

project_dir = Path(".")
download_workdir = project_dir / lcc.download_workdir_name
construct_counts_workdir = project_dir / lcc.construct_counts_workdir_name
qc_workdir = project_dir / lcc.qc_workdir_name
transform_workdir = project_dir / lcc.transform_workdir_name
type_cells_workdir = project_dir / lcc.type_cells_workdir_name
harmonize_metadata_workdir = project_dir / lcc.harmonize_metadata_workdir_name

def find_workdir_anndata(workdir: Path, name: str):
    matches = [p for p in workdir.iterdir() if p.is_file() and p.name == name]
    if len(matches) != 1:
        raise ValueError(f"Unable to find suitable anndata in {workdir}")
    return matches[0]

@main.command("download")
@click.argument("action", type=click.Choice(stepwise_actions))
@click.option("--gse-id", type=str)
def download(action: list[StepwiseAction], gse_id: str):
    if action == "run":

        subseries_ids = get_subseries_ids(gse_id)
        if len(subseries_ids) != 0:
            click.secho(
                    f"Detected SuperSeries   → {gse_id} contains {subseries_ids}"
                    "Proceed by constructing a directory and running the tool on each individual accession.",
                fg="red")
            return

        srp_cache: dict[str, pd.DataFrame] = {}
        click.secho(f"[download/run] Processing {gse_id}", bold=True)

        download_workdir.mkdir(exist_ok=True)
        construct_study_metadata(
            gse_id,
            download_workdir / lcc.metadata_file_name,
            srp_df_cache=srp_cache,
        )
        click.secho("  ↳ metadata ok", fg="green")

        download_gse_supps(
            gse_id,
            download_workdir / lcc.supp_data_dir_name,
        )
        click.secho("  ↳ supp files ok", fg="green")

        (download_workdir / lcc.external_id_file_name).write_text(gse_id)

        click.echo(
            f"[download/run] IMPORTANT: paste {gse_id} paper text "
            f"into {download_workdir / lcc.paper_text_file_name}"
        )
        click.echo(
            f"[download/run] IMPORTANT: paste {gse_id} paper URL "
            f"into {download_workdir / lcc.paper_url_file_name}"
        )

    elif action == "validate":
        for n in {lcc.metadata_file_name, lcc.paper_text_file_name, lcc.supp_data_dir_name}: 
            assert (project_dir / n).exists()
    else:
        raise ValueError(f"Invalid value {action}. Choose from {stepwise_actions}")

@main.command("construct-counts")
@click.argument("action", type=click.Choice(stepwise_actions))
def construct_counts(action: list[StepwiseAction]):

    if action == "run":
        supp_data_dir = download_workdir / lcc.supp_data_dir_name
        metadata_file = download_workdir / lcc.metadata_file_name
        paper_text_file = download_workdir / lcc.paper_text_file_name
        for n in {supp_data_dir, metadata_file, paper_text_file}:
            assert n.exists(), f"{n.name} does not exist"

        print("[construct-counts/run] Starting count matrix construction")
        _construct_counts(
            supp_data_dir,
            paper_text_file,
            metadata_file,
            construct_counts_workdir
        )

    elif action == "validate":
        assert (construct_counts_workdir / lcc.construct_counts_adata_name).exists()
    else:
        raise ValueError(f"Invalid value {action}. Choose from {stepwise_actions}")

@main.command("qc")
@click.argument("action", type=click.Choice(stepwise_actions))
def qc(action: list[StepwiseAction]):

    if action == "run":
        anndata_file = find_workdir_anndata(construct_counts_workdir,
                                            lcc.construct_counts_adata_name)
        print("[qc/run] Reading AnnData")
        adata = sc.read_h5ad(anndata_file)

        metadata_file = download_workdir / lcc.metadata_file_name
        paper_text_file = download_workdir / lcc.paper_text_file_name
        for n in {anndata_file, metadata_file, paper_text_file}:
            assert n.exists(), f"{n.name} does not exist"

        print("[qc/run] Starting count matrix construction")
        qc_and_filter(adata, paper_text_file, metadata_file, qc_workdir)

    elif action == "validate":
        assert (qc_workdir / lcc.qc_adata_name).exists()
    else:
        raise ValueError(f"Invalid value {action}. Choose from {stepwise_actions}")


@main.command("transform")
@click.argument("action", type=click.Choice(stepwise_actions))
def transform(action: list[StepwiseAction]):

    if action == "run":
        anndata_file = find_workdir_anndata(qc_workdir, lcc.qc_adata_name)

        print("[transform/run] Reading AnnData")
        adata = sc.read_h5ad(anndata_file)
        print("[transform/run] Starting counts transformation")
        transform_counts(adata, transform_workdir, use_scrublet=False)

    elif action == "validate":
        assert (transform_workdir / lcc.transform_adata_name).exists()
    else:
        raise ValueError(f"Invalid value {action}. Choose from {stepwise_actions}")

@main.command("type-cells")
@click.argument("action", type=click.Choice(stepwise_actions))
def type_cells(action: list[StepwiseAction]):

    if action == "run":
        anndata_file = find_workdir_anndata(transform_workdir, lcc.transform_adata_name)

        print("[type-cells/run] Reading AnnData")
        adata = sc.read_h5ad(anndata_file)
        print("[type-cells/run] Starting cell typing workflow")
        _type_cells(adata, type_cells_workdir)

    elif action == "validate":
        assert (type_cells_workdir / lcc.type_cells_adata_name).exists()
    else:
        raise ValueError(f"Invalid value {action}. Choose from {stepwise_actions}")

@main.command("harmonize-metadata")
@click.argument("action", type=click.Choice(stepwise_actions))
def harmonize_metadata(action: list[StepwiseAction]):

    if action == "run":
        anndata_file = find_workdir_anndata(type_cells_workdir, lcc.type_cells_adata_name)

        metadata_file = download_workdir / lcc.metadata_file_name
        paper_text_file = download_workdir / lcc.paper_text_file_name
        for n in {metadata_file, paper_text_file}:
            assert n.exists(), f"{n.name} does not exist"

        print("[harmonize-metadata/run] Reading AnnData")
        adata = sc.read_h5ad(anndata_file)
        print("[harmonize-metadata/run] Starting metadata harmonization workflow")
        _harmonize_metadata(adata, metadata_file, paper_text_file, harmonize_metadata_workdir)

    elif action == "validate":
        assert (harmonize_metadata_workdir / lcc.harmonize_metadata_adata_name).exists()
    else:
        raise ValueError(f"Invalid value {action}. Choose from {stepwise_actions}")

@main.command("publish-build")
@click.option("--anndata-path", type=click.Path(exists=True, path_type=Path))
def publish_build(anndata_path: Path):

        paper_text_file = Path(download_workdir / lcc.paper_text_file_name)
        paper_url_file = Path(download_workdir / lcc.paper_url_file_name)
        external_id_file = Path(download_workdir / lcc.external_id_file_name)

        assert paper_text_file.exists()
        assert paper_url_file.exists()
        assert external_id_file.exists()

        adata = sc.read_h5ad(anndata_path)
        is_error, tags = lint_anndata(adata)

        build_publish_data(
            paper_text_file.read_text(),
            paper_url_file.read_text(),
            external_id_file.read_text(),
            adata,
            Path(lcc.publish_workdir_name),
            tags
        )


@main.command("publish-upload")
@click.option("--anndata-path", type=click.Path(exists=True, path_type=Path))
@click.option("--latch-dest", type=str)
def publish_upload(anndata_path: Path, latch_dest: str):

        full_latch_dest = f'{latch_dest}/{anndata_path.name}'
        # latch_cp([str(anndata_path.resolve())], full_latch_dest, progress=Progress.tasks, verbose=False, expand_globs=False)

        res = get_node_data(full_latch_dest)
        node_id = res.data[full_latch_dest].id
        print(f'retrieved node id {node_id} for {full_latch_dest}')

        build_publish_queries(
            node_id,
            Path(lcc.publish_workdir_name) / lcc.publish_build_info_file_name,
            Path(lcc.publish_workdir_name) / "queries.sql",
        )

    # upload portal
    # send email


@main.command("convert")
def convert_scanpy_to_seurat(action: list[StepwiseAction]):
    ...

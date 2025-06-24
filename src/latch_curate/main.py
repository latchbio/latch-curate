from pathlib import Path
from typing import Literal
from textwrap import dedent

import scanpy as sc
import click
import pandas as pd
import json

from latch.ldata._transfer.node import get_node_data
from latch.ldata._transfer.progress import Progress
from latch_cli.services.cp.main import cp as latch_cp

from latch_curate.constants import latch_curate_constants as lcc
from latch_curate.download import construct_study_metadata, download_gse_supps, get_subseries_ids
from latch_curate.construct import construct_counts as _construct_counts
from latch_curate.qc import qc_and_filter
from latch_curate.transform import transform_counts
from latch_curate.cell_typing import type_cells as _type_cells
from latch_curate.harmonize import harmonize_metadata as _harmonize_metadata
from latch_curate.publish.build import build_publish_report, build_publish_data
from latch_curate.lint.linter import lint_anndata
from latch_curate.publish.queries import build_publish_queries
from latch_curate.publish.email_utils import EmailRecipient, send_email_to_authors
from latch_curate.utils import write_anndata
from latch_curate.llm_utils import get_token_count

StepwiseAction = Literal["run", "review"]
stepwise_actions: list[StepwiseAction] = ["run", "review"]

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
publish_workdir = project_dir / lcc.publish_workdir_name

def find_workdir_anndata(workdir: Path, name: str):
    matches = [p for p in workdir.iterdir() if p.is_file() and p.name == name]
    if len(matches) != 1:
        raise ValueError(f"Unable to find suitable anndata in {workdir}")
    return matches[0]


@main.group("download")
def download():
    ...

@download.command("run")
@click.option("--gse-id", type=str)
def download_run(gse_id: str):
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

def check_download_files_exist() -> (Path, Path, Path):
    supp_data_dir = download_workdir / lcc.supp_data_dir_name
    paper_text_file = download_workdir / lcc.paper_text_file_name
    metadata_file = download_workdir / lcc.metadata_file_name
    for n in {supp_data_dir, metadata_file, paper_text_file}:
        assert n.exists(), f"{n.name} does not exist"
    return supp_data_dir, paper_text_file, metadata_file


@download.command("review")
def download_review():
    _, paper_text_path, metadata_file_path = check_download_files_exist()
    paths = [paper_text_path, metadata_file_path]

    file_counts = {}
    total_file_tokens = 0
    for p in paths:
        count = get_token_count(p.read_text())
        file_counts[p.name] = count
        total_file_tokens += count

    # todo(kenny): estimate for o4-mini
    overhead_tokens = 600

    total_tokens = total_file_tokens + overhead_tokens

    threshold = 200_000
    click.echo(f"Threshold: {threshold} tokens\n")

    click.echo("File token counts:")
    for name, cnt in file_counts.items():
        click.echo(f"  • {name:<20} {cnt:>6,} tokens")
    click.echo(f"\nEstimated prompt/system overhead: {overhead_tokens:,} tokens")
    click.echo(f"\n>>> TOTAL TOKENS: {total_tokens:,}\n")

    if total_tokens > threshold:
        raise click.ClickException(
            f"Total ({total_tokens:,}) exceeds limit ({threshold:,}).\n"
            "    • Consider truncating your input files,\n"
            "    • Removing repetitive or uninformative sections (eg.  citations)\n"
        )

    click.secho("All good — total tokens within limit.", fg="green")

@main.group("construct-counts")
def construct_counts():
    ...


@construct_counts.command(name="run")
def construct_counts_run():
    supp_data_dir, paper_text_file, metadata_file = check_download_files_exist()
    print("[construct-counts/run] Starting count matrix construction")
    _construct_counts(
        supp_data_dir,
        paper_text_file,
        metadata_file,
        construct_counts_workdir
    )
    assert (construct_counts_workdir / lcc.construct_counts_adata_name).exists()

# def construct_counts_review(query: str):
#     _, paper_text_file, metadata_file = check_download_files_exist()
#     print("[construct-counts/review] Starting review of construction context")
#     review_counts(
#         paper_text_file,
#         metadata_file,
#         construct_counts_workdir,
#         query,
#     )
#     return

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
def publish_build():

        paper_text_file = Path(download_workdir / lcc.paper_text_file_name)
        paper_url_file = Path(download_workdir / lcc.paper_url_file_name)
        external_id_file = Path(download_workdir / lcc.external_id_file_name)

        assert paper_text_file.exists()
        assert paper_url_file.exists()
        assert external_id_file.exists()

        adata = sc.read_h5ad(harmonize_metadata_workdir /
                            lcc.harmonize_metadata_adata_name)
        is_error, tags = lint_anndata(adata)

        build_publish_data(
           paper_text_file.read_text(),
           paper_url_file.read_text(),
           external_id_file.read_text(),
           adata,
           publish_workdir,
           tags
        )

        report_map = {
            "construct_counts": construct_counts_workdir / lcc.construct_counts_report_name,
            "qc": qc_workdir / lcc.qc_report_name,
            "transform": transform_workdir / lcc.transform_report_name,
            "type_cells": type_cells_workdir / lcc.type_cells_report_name,
            "harmonize_metadata": harmonize_metadata_workdir / lcc.harmonize_metadata_report_name,
        }

        build_publish_report(
            publish_workdir,
            report_map
        )
        write_anndata(adata, publish_workdir, lcc.publish_adata_name)


@main.command("publish-upload")
@click.option("--latch-dest", type=str)
def publish_upload(latch_dest: str):

        publish_data = json.loads((publish_workdir /
                                   lcc.publish_build_info_file_name).read_text())

        display_name = publish_data['info']['display_name']
        full_latch_dest = f'{latch_dest}/{display_name}'

        latch_cp([str(publish_workdir.resolve())], full_latch_dest, progress=Progress.tasks, verbose=False, expand_globs=False)

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

@main.command("publish-email")
def publish_email():

    publish_data = json.loads((publish_workdir /
                               lcc.publish_build_info_file_name).read_text())

    dataset_info = publish_data['info']

    corr_author_names = dataset_info['corresponding_author_names']
    corr_author_emails = dataset_info['corresponding_author_emails']
    paper_title = dataset_info['paper_title']
    display_name = dataset_info['display_name']
    cell_count = dataset_info['cell_count']

    recipients = [EmailRecipient(email=email,name=name) for email, name in
                  zip(corr_author_emails, corr_author_names)]


    subject = f"Automated LLM Curation of {paper_title}"
    greeting = ""

    assert len(recipients) > 0

    recipient_names = [x.name for x in recipients]
    if len(recipient_names) == 1:
        greeting = f"Dear {recipient_names[0]},"
    elif len(recipient_names) == 2:
        greeting = f"Dear {recipient_names[0]} and {recipient_names[1]},"
    else:
        all_but_last = ", ".join(recipient_names[:-1])
        last = recipient_names[-1]
        greeting = f"Dear {all_but_last}, and {last},"

    body = dedent(
        f"""\
<html>
  <body>
    <p>{greeting}</p>

    <p>The {cell_count} single cells associated with <strong>{paper_title}</strong> have been
    automatically curated with an agential LLM framework and uploaded to
    <a href="https://console.latch.bio/datasets" target="_blank">LatchBio</a>.<p>

    <p>A report outlining e.g. quality control, cell typing,
    study metadata is attached to this email. It also provides an overview
    of the curation mechanics and how our team was able to automate this
    process.</p>

    <p>The dataset is currently being sold on Latch as {display_name.strip()}. 
    Companies in industry will use it to create large cell atlases and train new
    foundation models - towards understanding, treating and eventually
    curing complex disease.</p>

    <p>Our intention is to redirect a majority of this revenue back to the
    authors to fund more basic research. Please respond to this email to
    participate in revenue sharing or for any other inquiries.</p>

    <p>Thank you for your contributions to science and the future of
    engineering biology.</p>

    <p>Kenny Workman<br/>
    CTO | LatchBio</p>

    <hr/>

    <p><small>This message and curation results are entirely automated. Excuse any errors.</small></p>
  </body>
</html>
"""
    )

    send_email_to_authors(subject, body, recipients, publish_workdir / lcc.publish_report_name)


@main.command("convert")
def convert_scanpy_to_seurat(action: list[StepwiseAction]):
    ...

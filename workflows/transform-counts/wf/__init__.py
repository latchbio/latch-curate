from pathlib import Path


import scanpy as sc

from latch_curate.transform import transform_counts

from latch.resources.workflow import workflow
from latch.types.directory import LatchOutputDir
from latch.types.file import LatchFile
from latch.types.metadata import LatchAuthor, LatchMetadata, LatchParameter
from latch.resources.tasks import medium_task, large_task

metadata = LatchMetadata(
    display_name="Transform Counts",
    author=LatchAuthor(
        name="Kenny Workman",
    ),
    parameters={
        "input_h5ad": LatchParameter(
            display_name="Input H5AD",
            batch_table_column=True,  # Show this parameter in batched mode.
        ),
        "use_scrublet": LatchParameter(
            display_name="Use Scrublet",
        ),
        "output_directory": LatchParameter(
            display_name="Output Directory",
            batch_table_column=True,  # Show this parameter in batched mode.
        ),
    },
)

large_metadata = LatchMetadata(
    display_name="Transform Counts (Big Computer)",
    author=LatchAuthor(
        name="Kenny Workman",
    ),
    parameters={
        "input_h5ad": LatchParameter(
            display_name="Input H5AD",
            batch_table_column=True,  # Show this parameter in batched mode.
        ),
        "output_directory": LatchParameter(
            display_name="Output Directory",
            batch_table_column=True,  # Show this parameter in batched mode.
        ),
    },
)


@medium_task
def transform_counts_task(input_h5ad: LatchFile, output_directory:
                          LatchOutputDir, use_scrublet: bool) -> (LatchFile, LatchFile):
    input_path = Path(input_h5ad.local_path)

    adata = sc.read_h5ad(input_path)
    adata, html_path = transform_counts(adata, use_scrublet=use_scrublet)

    output_h5ad_name = f"{input_path.stem}.transformed.h5ad"
    adata.write_h5ad(output_h5ad_name)

    h5ad_latch = LatchFile(
        str(Path(output_h5ad_name).resolve()),
        f"{output_directory.remote_directory}/{output_h5ad_name}"
    )
    html_fname = Path(html_path).name
    html_latch = LatchFile(
        str(Path(html_path).resolve()),
        f"{output_directory.remote_directory}/{html_fname}"
    )

    return h5ad_latch, html_latch

@large_task
def transform_counts_task_large(input_h5ad: LatchFile, output_directory:
                                LatchOutputDir, use_scrublet: bool = True) -> (LatchFile, LatchFile):
    input_path = Path(input_h5ad.local_path)

    adata = sc.read_h5ad(input_path)
    adata, html_path = transform_counts(adata, use_scrublet)

    output_h5ad_name = f"{input_path.stem}.transformed.h5ad"
    adata.write_h5ad(output_h5ad_name)

    h5ad_latch = LatchFile(
        str(Path(output_h5ad_name).resolve()),
        f"{output_directory.remote_directory}/{output_h5ad_name}"
    )
    html_fname = Path(html_path).name
    html_latch = LatchFile(
        str(Path(html_path).resolve()),
        f"{output_directory.remote_directory}/{html_fname}"
    )

    return h5ad_latch, html_latch

@workflow(metadata)
def transform_counts_workflow(
        input_h5ad: LatchFile, output_directory: LatchOutputDir, use_scrublet:
        bool = True
) -> (LatchFile, LatchFile):
    return transform_counts_task(input_h5ad=input_h5ad,
                                 output_directory=output_directory,
                                 use_scrublet=use_scrublet)

# @workflow(large_metadata)
# def transform_counts_large(
#     input_h5ad: LatchFile, output_directory: LatchOutputDir
# ) -> (LatchFile, LatchFile):
#     return transform_counts_task_large(input_h5ad=input_h5ad, output_directory=output_directory)

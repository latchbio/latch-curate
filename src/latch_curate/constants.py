from dataclasses import dataclass

@dataclass(frozen=True)
class LatchCurateConstants:
    pkg_name: str = "latch-curate"

    pkg_version_cache_path: str = "latch-curate/cached-version.txt"
    openai_api_key_path: str = "latch-curate/openai_api_key.txt"

    metadata_file_name: str = "study_metadata.txt"
    paper_text_file_name: str = "paper_text.txt"
    supp_data_dir_name: str = "supp_data"

    construct_counts_adata_name: str = "counts.h5ad"
    construct_counts_workdir_name: str = "construct_counts"

    qc_adata_name: str = "qc.h5ad"
    qc_workdir_name: str = "qc"
    qc_report_name: str = "qc_report.html"

    transform_adata_name: str = "transform.h5ad"
    transform_workdir_name: str = "transform"
    transform_report_name: str = "transform.html"

    type_cells_adata_name: str = "type_cells.h5ad"
    type_cells_workdir_name: str = "type_cells"
    type_cells_report_name: str = "type_cells.html"

    harmonize_metadata_adata_name: str = "harmonize_metadata.h5ad"
    harmonize_metadata_workdir_name: str = "harmonize_metadata"
    harmonize_metadata_report_name: str = "harmonize_metadata.html"


latch_curate_constants = LatchCurateConstants()

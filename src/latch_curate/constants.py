from dataclasses import dataclass

@dataclass(frozen=True)
class LatchCurateConstants:
    pkg_name: str = "latch-curate"
    nucleus_url: str = "https://nucleus.latch.bio/infer"
    # nucleus_url: str = "http://localhost:5000"

    get_cell_count_endpoint = "get-cell-count"
    get_fixed_qc_thresholds_endpoint = "get-fixed-qc-thresholds"
    get_cell_types_endpoint = "get-cell-types"
    get_harmonized_metadata_endpoint = "get-harmonized-metadata"

    pkg_version_cache_path: str = "latch-curate/cached-version.txt"
    openai_api_key_path: str = "latch-curate/openai_api_key.txt"

    external_id_file_name: str = "external_id.txt"
    metadata_file_name: str = "study_metadata.txt"
    paper_text_file_name: str = "paper_text.txt"
    paper_url_file_name: str = "paper_url.txt"
    supp_data_dir_name: str = "supp_data"

    download_workdir_name: str = "download"

    construct_counts_adata_name: str = "counts.h5ad"
    construct_counts_workdir_name: str = "construct_counts"
    construct_counts_report_name: str = "construct_counts.html"

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

    publish_workdir_name: str = "publish"
    publish_adata_name: str = "publish.h5ad"
    publish_report_name: str = "publish.html"

    publish_build_info_file_name: str = "build.json"


latch_curate_constants = LatchCurateConstants()

from pathlib import Path
from latch_curate.download.metadata import construct_study_metadata
from latch_curate.download.supps import download_gse_supps
from latch_curate.construct.harness import construct_counts
from latch_curate.cell_typing import type_cells
from latch_curate.qc import qc_and_filter
from latch_curate.transform import transform_counts


gse_id = "GSE252545"

counts_path = Path("counts")
data_path = counts_path / "data"
data_path.mkdir(parents=True, exist_ok=True)

# 1/ Download paper text + supps
# - use browser automation to avoid copying and paseting
# - explore PDFs

construct_study_metadata(gse_id)
download_gse_supps(gse_id, data_path)

# 2/ Construct counts using paper text, study metadata + various tools

adata_file = construct_counts(
        data_path, 
        Path('paper_text.txt'),
        Path('study_metadata.txt'),
        Path('counts'), 
)

# quit()

adata = adata_file()

# 3/ QC with fixed thresholds inferred from paper + metadata and adaptive
# thresholds built from quantile tables

qc_report_path = qc_and_filter(adata, Path('study_metadata.txt'), Path('paper_text.txt'))

# 4/ Count transformation. Raw counts to UMAP. No LLM (yet).

transform_report_path = transform_counts(adata, use_scrublet=False)

# 5/ Dump diff-exp genes into JSON. Combine with paper + metadata to type
# clusters.

cell_typing_report_path = type_cells(adata)

# 6/ Metadata harmonization using all paper text and study metadata against
# whatever controlled vocabulary.

# TODO

# 7/ Linting + Seurat conversion

# Latch workflows. TODO: allow awaitable workflow execution from python

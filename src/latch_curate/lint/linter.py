from anndata import AnnData

from latch_curate.lint.vocab.cl import validate_cl
from latch_curate.lint.vocab.mondo import validate_mondo
from latch_curate.lint.vocab.sample_site import validate_sample_site
from latch_curate.lint.vocab.sequencing_platform import validate_efo
from latch_curate.lint.vocab.uberon import validate_uberon
from latch_curate.lint.vocab.common import unknown_val
from latch_curate.publish.build import Tag

obs_names = [
    "latch_cell_type_lvl_1",
    "latch_disease",
    "latch_tissue",
    "latch_sample_id",
    "latch_subject_id",
    "latch_condition",
    "latch_sample_site",
    "latch_sequencing_platform",
    "latch_organism",
]

controlled_obs_map = {
    "latch_cell_type_lvl_1": validate_cl,
    "latch_disease": validate_mondo,
    "latch_tissue": validate_uberon,
    "latch_sample_site": validate_sample_site,
    "latch_sequencing_platform": validate_efo,
}

obs_key_to_orion_metadata = {
    "latch_cell_type_lvl_1": "cell_type",
    "latch_disease": "disease",
    "latch_tissue": "tissue",
    "latch_sequencing_platform": "assay",
}

var_names = [
    "gene_symbols",
    "highly_variable",
    "means",
    "dispersions",
    "dispersions_norm",
    "highly_variable_nbatches",
    "highly_variable_intersection",
]

uns_names = [
    "hvg",
    "pca",
    "latch_sample_id_colors",
    "neighbors",
    "umap",
    "leiden_res_0.50",
]

obsm_names = ["X_pca", "X_umap"]

def lint_anndata(adata: AnnData) -> (bool, list[Tag]):
    tags = []
    is_error = False
    for obs_name in obs_names:
        if obs_name not in adata.obs.columns:
            is_error = True
            print(f">>> [obs] {obs_name} absent.")

    for var_name in var_names:
        if var_name not in adata.var.columns:
            is_error = True
            print(f">>> [var] {var_name} absent.")

    uns_keys = set(adata.uns.keys())
    for uns_name in uns_names:
        if uns_name not in uns_keys:
            is_error = True
            print(f">>> [uns] {uns_name} absent.")

    for obsm_name in obsm_names:
        if obsm_name not in adata.obsm:
            is_error = True
            print(f">>> [obsm] {obsm_name} absent.")

    for k, validater in controlled_obs_map.items():
        for x in set(adata.obs[k]):
            is_valid = validater(x)
            if not is_valid:
                is_error = True
                print(f">>> [{k}] {x} not valid.")
            else:
                if k in obs_key_to_orion_metadata:
                    if x == unknown_val:
                        continue
                    name, obo_id = x.split("/")
                    tags.append(Tag(metadata_type=obs_key_to_orion_metadata[k], name=name, obo_id=obo_id))
                print(f"[{k}] {x} valid")

    return is_error, tags

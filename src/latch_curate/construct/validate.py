import re

import numpy as np
import scipy.sparse as sp
import anndata as ad

def record_and_assert(validation_log: list[tuple[str, str]], cond: bool, msg: str):
    validation_log.append(("pass" if cond else "fail", msg))
    assert cond, msg

ensembl_pattern = re.compile(r"ENS[A-Z0-9]{0,5}G\d{11}(?:\.\d+)?$")

def validate_counts_object(adata: ad.AnnData) -> list[tuple[str, str]]:

    validation_log: list[tuple[str, str]] = []

    record_and_assert(validation_log, 'latch_sample_id' in adata.obs, "obs contains latch_sample_id")

    record_and_assert(validation_log, all(map(bool, map(ensembl_pattern.match, adata.var_names))),
            "var index are Ensembl IDs")

    record_and_assert(validation_log, 'gene_symbols' in adata.var.columns, "var contains gene_symbols")
    record_and_assert(validation_log, adata.var['gene_symbols'].is_unique, "gene_symbols unique")

    record_and_assert(validation_log, adata.obs_names.is_unique, "obs index unique")
    record_and_assert(validation_log, adata.var_names.is_unique, "var index unique")

    X = adata.X
    vals = X.data if sp.issparse(X) else np.asarray(X)
    record_and_assert(validation_log, (vals >= 0).all(), "counts non-negative")
    record_and_assert(validation_log, np.allclose(vals, vals.round()), "counts are integers")

    for col in adata.obs.columns:
        if col.startswith("author_"):
            series = adata.obs[col]
            record_and_assert(validation_log, ~series.isna().all(), f"{col} not all NaN")

    # n_obs = adata.n_obs
    # lower = target_cell_count - tol
    # upper = target_cell_count + tol
    # cond = (lower <= n_obs <= upper)
    # record_and_assert(
    #     validation_log,
    #     cond,
    #     f"n_obs={n_obs} within ±{tol} of target {target_cell_count}"
    # )

    return validation_log

# scrna_utils.py

import warnings
from pathlib import Path
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
import json
import re

# todo(kenny): pull out
ENSEMBL_REGEX = re.compile(r"ENS[A-Z0-9]{0,5}G\d{11}(?:\.\d+)?$")


# suppress pandas dtype warnings
warnings.simplefilter("ignore", pd.errors.DtypeWarning)

ensembl_map: dict[str, str] = {}
with open('ensembl_map.json') as f:
    ensembl_map = json.load(f)

def is_valid_ensembl(id_: str) -> bool:
    """True if `id_` matches the Ensembl-gene pattern enforced by the validator."""
    return ENSEMBL_REGEX.fullmatch(id_) is not None

def ensembl_to_symbol(id_: str, mapping: dict[str, str] = ensembl_map) -> str | None:
    """Return the gene symbol for a valid Ensembl ID, or None if not found/invalid."""
    if not is_valid_ensembl(id_):
        return None
    if "_reverse_map" not in ensembl_to_symbol.__dict__:
        ensembl_to_symbol._reverse_map = {v: k for k, v in mapping.items()}
    return ensembl_to_symbol._reverse_map.get(id_)


def symbol_to_ensembl(sym: str, mapping: dict[str, str] = ensembl_map) -> str | None:
    """Return the Ensembl ID for a gene symbol, or None if unknown/invalid."""
    return mapping.get(sym.upper())


def convert_and_swap_symbol_index(
    adata: AnnData,
    mapping: dict[str, str],
) -> AnnData:
    """
    Replace var.index gene symbols with Ensembl IDs based on `mapping`.
    Keeps only genes that map, stores original symbols in var['gene_symbols'],
    and removes duplicate var_names after mapping.
    """
    symbols = pd.Index(adata.var.index)
    ensembl_ids = symbols.to_series().map(mapping)
    mask = ~ensembl_ids.isna()

    # filter to mapped genes
    adata = adata[:, mask].copy()

    # record original symbols and swap index
    adata.var["gene_symbols"] = adata.var.index
    adata.var.index = ensembl_ids[mask]  # now Ensembl IDs

    # remove any duplicate Ensembl IDs
    unique_mask = ~adata.var_names.duplicated()
    adata = adata[:, unique_mask].copy()

    return adata


def ingest_csv_data(
    counts_path: Path,
    annotation_path: Path,
    mapping: dict[str, str],
    delimiter: str = "\t",
    metadata_keys: set[str] = {"sample_stuff"},
) -> AnnData:
    """
    Read a counts CSV/TSV and a metadata CSV/TSV, build an AnnData.

    - counts_path: TSV with barcodes in columns (after first), genes in rows.
    - annotation_path: TSV with cell metadata; must have 'cellname' column.
    - mapping: symbol→Ensembl mapping passed to convert_and_swap_symbol_index.
    - delimiter: field delimiter for both files (default tab).
    - metadata_keys: which columns in annotation to pull into obs (prefix 'author_').
    """
    # --- counts ---
    df = pd.read_csv(counts_path, delimiter=delimiter, engine="c", memory_map=True)
    cell_names = df.columns[1:].tolist()

    # drop header row if present, extract gene names and expression
    gene_df = df.iloc[1:, :]
    gene_names = gene_df.iloc[:, 0].tolist()
    expr_data = gene_df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")

    # build expression matrix
    X = expr_data.T
    X.index = cell_names
    X.columns = gene_names

    # --- metadata ---
    meta = pd.read_csv(annotation_path, delimiter=delimiter, engine="c", memory_map=True)
    meta = meta.set_index("cellname").reindex(cell_names)

    obs_dict = {
        f"author_{k}": meta[k]
        for k in metadata_keys
        if k in meta.columns
    }
    obs = pd.DataFrame(obs_dict, index=cell_names)

    # assemble AnnData
    adata = AnnData(X=X, obs=obs)
    adata = convert_and_swap_symbol_index(adata, mapping)
    return adata


def reindex_and_fill_with_zeros(
    adata: AnnData,
    new_features: pd.Index,
    symbol_map: dict[str, str],
) -> AnnData:
    """
    Reindex `adata` to the given `new_features`, filling missing genes with zeros.
    Updates var to include 'gene_symbols' from `symbol_map`.
    """
    # densify if needed
    X_dense = adata.X.todense() if sp.issparse(adata.X) else adata.X
    df = pd.DataFrame(X_dense, index=adata.obs_names, columns=adata.var_names)

    # reindex columns → new union of features
    df_reindexed = df.reindex(columns=new_features, fill_value=0)

    # rebuild var
    new_var = adata.var.reindex(new_features).copy()
    new_var["gene_symbols"] = new_var.index.to_series().map(symbol_map)

    # build new AnnData
    return AnnData(
        X=sp.csr_matrix(df_reindexed.values),
        obs=adata.obs.copy(),
        var=new_var,
    )


def reindex_and_fill_list(adatas: list[AnnData]) -> list[AnnData]:
    """
    Given a list of AnnData objects, compute the union of all feature indices,
    then reindex each to that union (filling zeros) and return the new list.
    """
    # union of all var indices
    union_features = pd.Index([]).union_many([ad.var.index for ad in adatas])
    print(f"Size of union index: {len(union_features)}")

    # build symbol_map from first occurrence in each adata
    symbol_map: dict[str, str] = {}
    for feature in union_features:
        for ad in adatas:
            if feature in ad.var.index:
                symbol_map[feature] = ad.var.loc[feature, "gene_symbols"]
                break

    # reindex each
    new_adatas = []
    for ad in adatas:
        new_ad = reindex_and_fill_with_zeros(ad, union_features, symbol_map)
        new_adatas.append(new_ad)
        print(">>> Reindexed + filled:", new_ad)
    return new_adatas

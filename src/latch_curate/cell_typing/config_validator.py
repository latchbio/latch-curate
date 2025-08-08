from pathlib import Path
from dataclasses import dataclass
import yaml
from anndata import AnnData


@dataclass
class CellTypeVocabularyEntry:
    name: str
    ontology_id: str


@dataclass
class CellTypingConfig:
    cell_type_column: str
    cluster_column: str
    vocabulary: list[CellTypeVocabularyEntry]
    marker_genes: dict[str, list[str]]


@dataclass
class CellTypingTag:
    cell_type: str
    ontology_id: str


def parse_cell_typing_config(path: Path) -> CellTypingConfig:
    with open(path) as f:
        data = yaml.safe_load(f)

    vocab_entries = []
    for entry in data.get('vocabulary', []):
        vocab_entries.append(CellTypeVocabularyEntry(
            name=entry['name'],
            ontology_id=entry['ontology_id']
        ))

    config = CellTypingConfig(
        cell_type_column=data.get('cell_type_column', 'latch_cell_type_lvl_1'),
        cluster_column=data.get('cluster_column', 'leiden_res_0.50'),
        vocabulary=vocab_entries,
        marker_genes=data.get('marker_genes', {})
    )

    return config


def validate_cell_typing_config(config: CellTypingConfig, adata: AnnData) -> tuple[bool, list[str]]:
    errors = []

    if config.cluster_column not in adata.obs.columns:
        errors.append(f"Cluster column '{config.cluster_column}' not found in adata.obs")

    if config.cell_type_column in adata.obs.columns:
        unique_types = adata.obs[config.cell_type_column].dropna().unique()
        valid_types = {v.name for v in config.vocabulary}
        valid_types_with_ontology = {f"{v.name}/{v.ontology_id}" for v in config.vocabulary}
        
        for cell_type in unique_types:
            cell_type_name = cell_type.split('/')[0] if '/' in cell_type else cell_type
            if cell_type not in valid_types and cell_type not in valid_types_with_ontology and cell_type_name not in valid_types:
                errors.append(f"Cell type '{cell_type}' not in configured vocabulary")

    if 'gene_symbols' in adata.var.columns:
        available_genes = set(adata.var['gene_symbols'])
        for cell_type, genes in config.marker_genes.items():
            missing_genes = [g for g in genes if g not in available_genes]
            if missing_genes:
                print(f"Warning: {len(missing_genes)} marker genes for '{cell_type}' not found in data")

    return len(errors) == 0, errors


def extract_cell_typing_tags(adata: AnnData, config: CellTypingConfig) -> list[CellTypingTag]:
    tags = []

    if config.cell_type_column not in adata.obs.columns:
        return tags

    vocab_map = {v.name: v.ontology_id for v in config.vocabulary}

    unique_types = adata.obs[config.cell_type_column].dropna().unique()

    for cell_type in unique_types:
        if cell_type in vocab_map:
            tags.append(CellTypingTag(
                cell_type=cell_type,
                ontology_id=vocab_map[cell_type]
            ))

    return tags

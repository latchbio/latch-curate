from pathlib import Path
from dataclasses import dataclass
from typing import Optional
import yaml

@dataclass
class Vocab:
    type: str
    name: Optional[str] = None
    values: Optional[list[str]] = None

@dataclass
class VariableDefinition:
    name: str
    description: str
    vocab: Vocab

def var_to_json(var: VariableDefinition) -> dict:
    vocab_data = {"type": var.vocab.type}
    if var.vocab.name is not None:
        vocab_data["name"] = var.vocab.name
    if var.vocab.values is not None:
        vocab_data["values"] = var.vocab.values

    return {
        "name": var.name,
        "description": var.description,
        "vocab": vocab_data,
    }


def parse_metadata_yaml(path: Path) -> list[VariableDefinition]:
    with open(path) as f:
        data = yaml.safe_load(f)
    vars_list = []
    for entry in data.get('variables', []):
        vocab_data = entry.get('vocab', {})
        vocab = Vocab(
            type=vocab_data.get('type'),
            name=vocab_data.get('name'),
            values=vocab_data.get('values')
        )
        var_def = VariableDefinition(
            name=entry.get('name'),
            description=entry.get('description'),
            vocab=vocab
        )
        vars_list.append(var_def)
    return vars_list

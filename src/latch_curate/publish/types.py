from dataclasses import dataclass

@dataclass
class Tag:
    metadata_type: str
    name: str
    obo_id: str

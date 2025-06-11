import re
from contextlib import contextmanager
from tempfile import TemporaryDirectory
from pathlib import Path

import GEOparse

@contextmanager
def _tmpdir():
    with TemporaryDirectory() as td:
        yield Path(td)

super_series_pattern = r"SuperSeries of: (GSE\d+)"

def get_subseries_ids(gse_id: str) -> list[str]:
    with _tmpdir() as tmp:
        gse = GEOparse.get_GEO(geo=gse_id, destdir=str(tmp))
    rels = gse.metadata.get("relation", [])
    pat = re.compile(super_series_pattern, flags=re.IGNORECASE)
    return [m.group(1) for r in rels for m in [pat.match(r)] if m]

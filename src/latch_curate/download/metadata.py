from pathlib import Path
from textwrap import dedent
import traceback
import sys

from pysradb.sraweb import SRAweb
import GEOparse
import pandas as pd

from latch_curate.utils import df_to_str
sradb = SRAweb()

def remove_dir_recursive(p: Path):
    if not p.exists():
        return

    for child in p.iterdir():
        if child.is_dir():
            remove_dir_recursive(child)
        else:
            child.unlink()
    p.rmdir()

def download_srp_metadata(gse_id: str, srp_df_cache: dict[str, pd.DataFrame]):
    try:
        df = sradb.gse_to_srp(
            [gse_id],
        )
    except Exception:
        return

    assert df.shape[0] == 1
    assert df.columns[1] == 'study_accession'

    srp_id = df['study_accession'][0]

    if srp_id in srp_df_cache:
        return srp_df_cache[srp_id]

    df = sradb.sra_metadata(
        srp_id,
        detailed=True
    )
    srp_df_cache[srp_id] = df.copy()
    return df

def download_gse_metadata(gse_id: str):
    tmp_dir = Path("tmp")
    try:
        gse = GEOparse.get_GEO(geo=gse_id, destdir=str(tmp_dir.resolve()))
    finally:
        remove_dir_recursive(tmp_dir)

    all_keys = set()
    for gsm in gse.gsms.values():
        all_keys.update(gsm.metadata.keys())
    all_keys = sorted(list(all_keys))

    data = []
    for gsm_name, gsm in gse.gsms.items():
        row_dict = {"GSM": gsm_name}
        for key in all_keys:
            if key in gsm.metadata:
                row_dict[key] = "; ".join(gsm.metadata[key])
            else:
                row_dict[key] = ""
        data.append(row_dict)
    return pd.DataFrame(data)

_original_exit = sys.exit
def _fake_exit(code=0):
    raise RuntimeError(f"Intercepted sys.exit({code})")

def construct_study_metadata(gse_id: str, metadata_file: Path, srp_df_cache: dict[str, pd.DataFrame]):

    try:
        print("Downloading SRP metadata.")
        sys.exit = _fake_exit
        srp_df = download_srp_metadata(gse_id, srp_df_cache)
    except Exception:
        print(">>> Exception downloading SRP metadata: ")
        traceback.print_exc()
        srp_df = pd.DataFrame()
        print("Moving on")
    finally:
        sys.exit = _original_exit
    try:
        print("Downloading GSE metadata.")
        gse_df = download_gse_metadata(gse_id)
    except Exception:
        print(">>> Exception downloading GSE metadata: ")
        traceback.print_exc()
        gse_df = pd.DataFrame()
        print("Moving on")

    metadata_text = dedent(f"""
    <srp_metadata>
    {df_to_str(srp_df) if srp_df is not None else ""}
    </srp_metadata>

    <gse_metadata>
    {df_to_str(gse_df)}
    </gse_metadata>
    """)

    metadata_file.write_text(metadata_text)

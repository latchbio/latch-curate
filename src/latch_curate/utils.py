import io
import base64
from datetime import datetime, timedelta
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from anndata import AnnData

from latch_curate.constants import latch_curate_constants
from latch_curate.tinyrequests import get
from latch_curate.config import user_config

def _fig_to_base64(obj) -> str:

    if isinstance(obj, (list, tuple)) and obj:
        obj = obj[0]

    if hasattr(obj, "savefig"):
        fig = obj
    elif hasattr(obj, "figure"):
        fig = obj.figure
    elif hasattr(obj, "fig"):
        fig = obj.fig
    else:
        raise TypeError(f"Cannot extract Figure from object of type {type(obj)}")

    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("ascii")


def df_to_str(df: pd.DataFrame) -> str:
    return df.to_string(max_rows=None, max_cols=None, line_width=None)


def get_local_package_version() -> str:
    try:
        from importlib import metadata
    except ImportError:
        import importlib_metadata as metadata
    return metadata.version(latch_curate_constants.pkg_name)

def request_latest_package_version() -> str:
    resp = get(f"https://pypi.org/pypi/{latch_curate_constants.pkg_name}/json")
    version =  resp.json()["info"]["version"]
    user_config.package_version_cache_location.write_text(f"{version} {datetime.now().isoformat()}")
    print(user_config.package_version_cache_location)
    return version

def get_latest_package_version() -> str:
    version = None
    try:
        version, timestamp = user_config.package_version_cache_location.read_text().split(" ")
        if datetime.now() > datetime.fromisoformat(timestamp) + timedelta(days=1):
            version = request_latest_package_version()
    except FileNotFoundError:
        version = request_latest_package_version()
    return version


def write_anndata(adata: AnnData, workdir: Path, anndata_name: str):
    adata_path = workdir / anndata_name
    adata.write_h5ad(adata_path)
    print(f"anndata written to {adata_path.resolve()}")

def write_html_report(html: str, workdir: Path, report_name: str):
    output_html_path = workdir / report_name
    output_html_path.write_text(html)
    print(f"report written to {output_html_path.resolve()}")

def _df_to_html(df: pd.DataFrame) -> str:
    return (
        df.to_html(classes="quant", border=0, na_rep="", float_format="{:.2f}".format)
        .replace("<th>", "<th style='text-align:center'>")
        .replace("<td>", "<td style='text-align:center'>")
    )

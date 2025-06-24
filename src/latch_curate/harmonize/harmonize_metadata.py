from pathlib import Path
from textwrap import dedent
from html import escape
import json
from enum import Enum

from anndata import AnnData

from latch_curate.utils import write_html_report, write_anndata
from latch_curate.constants import latch_curate_constants as lcc
from latch_curate.config import user_config
from latch_curate.tinyrequests import post

class ControlledMetadataKeysEnum(str, Enum):
    latch_subject_id           = "latch_subject_id"
    latch_condition            = "latch_condition"
    latch_disease              = "latch_disease"
    latch_tissue               = "latch_tissue"
    latch_sample_site          = "latch_sample_site"
    latch_sequencing_platform  = "latch_sequencing_platform"
    latch_organism             = "latch_organism"

ControlledMetadataKeys = ControlledMetadataKeysEnum

def build_metadata_report(responses: dict[ControlledMetadataKeys, dict]) -> str:

    title = "Metadata Harmonization Report"

    rows_html: list[str] = []

    for key, payload in responses.items():
        ann = payload.get("annotations", {})
        reasoning_md = payload.get("reasoning", "")

        table_rows = "".join(
            f"<tr><td>{escape(sample_id)}</td><td>{escape(label)}</td></tr>"
            for sample_id, label in ann.items()
        )
        table_html = (
            "<table>"
            "<tr><th>Sample ID</th><th>Annotation</th></tr>"
            f"{table_rows}"
            "</table>"
        )

        reasoning_html = (
            "<h3>Chain&nbsp;of&nbsp;thought</h3>"
            f"<pre>{escape(reasoning_md)}</pre>"
        )

        rows_html.append(
            f"<h2>{escape(key)}</h2>"
            f"{table_html}"
            f"{reasoning_html}"
        )

    body = "\n".join(rows_html)

    return dedent(f"""\
        <!DOCTYPE html>
        <html lang="en"><head>
          <meta charset="utf-8">
          <title>{escape(title)}</title>
          <style>
            body {{ font-family: sans-serif; max-width: 900px; margin: auto; }}
            h2 {{ border-bottom: 1px solid #ccc; margin-top: 2em; }}
            table {{ border-collapse: collapse; width: 100%; margin-bottom: 1em; }}
            th, td {{ border: 1px solid #ddd; padding: 4px 8px; text-align: left; }}
            pre {{ white-space: pre-wrap; background: #fafafa; padding: 8px; 
                   border: 1px solid #eee; overflow-x: auto; }}
          </style>
        </head><body>
          <h1>{escape(title)}</h1>
          {body}
        </body></html>""")


def harmonize_metadata(
    adata: AnnData,
    metadata_file: Path,
    paper_text_file: Path,
    workdir: Path,
):
    workdir.mkdir(exist_ok=True)

    study_metadata = metadata_file.read_text()
    paper_text = paper_text_file.read_text()
    sample_list = list(set(adata.obs['latch_sample_id']))
    print(f"Harmonizing against sample_list from obs['latch_sample_id']: {sample_list}")

    response_dict = {}
    for k in [e.value for e in ControlledMetadataKeysEnum]:
        print(f"Requesting harmonized metadata from model for {k}")
        resp = post(
            f"{lcc.nucleus_url}/{lcc.get_harmonized_metadata_endpoint}",
            {
                "study_metadata": study_metadata,
                "paper_text": paper_text,
                "sample_list": json.dumps(sample_list),
                "metadata_key": k,
                "metadata": json.dumps({"step": "harmonize-metadata",
                                        "project": workdir.resolve().parent.name}),
                "session_id": -1
            },
            headers = {"Authorization": f"Latch-SDK-Token {user_config.token}"}
        )
        try:
            print(resp.json())
            data = resp.json()['data']
            annotations = data['annotations']
            reasoning = data['reasoning']
        except KeyError:
            raise ValueError(f'Malformed response data: {resp.json()}')

        response_dict[k] = {"annotations": annotations, "reasoning": reasoning}

    for k, v in response_dict.items():
        try:
            adata.obs[k] = adata.obs["latch_sample_id"].map(v["annotations"])
        except Exception as e:
            raise ValueError(f"Issue mapping {k}; {e.with_traceback()}")

    html = build_metadata_report(response_dict)
    write_html_report(html, workdir, lcc.harmonize_metadata_report_name)
    write_anndata(adata, workdir, lcc.harmonize_metadata_adata_name)

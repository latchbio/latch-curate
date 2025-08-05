from pathlib import Path
from textwrap import dedent
from html import escape
import yaml

from anndata import AnnData

from latch_curate.utils import write_html_report, write_anndata
from latch_curate.constants import latch_curate_constants as lcc
from latch_curate.config import user_config
from latch_curate.harmonize.schema import parse_metadata_yaml, var_to_json
from latch_curate.tinyrequests import post

annotation_dict_type = dict[str, [dict[str, str], str]]

def build_metadata_report(responses: annotation_dict_type) -> str:

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
    use_metadata: bool,
):
    workdir.mkdir(exist_ok=True)
    cache_path = workdir / lcc.harmonize_metadata_metadata_name

    annotation_dict: annotation_dict_type = {}
    if use_metadata:
        assert cache_path.exists(), f"Missing metadata file: {cache_path}"
        with cache_path.open() as f:
            annotation_dict = yaml.safe_load(f)
    else:
        try:
            var_defs = parse_metadata_yaml(user_config.metadata_schema_path)
        except Exception:
            raise ValueError('Malformed variable schema file.')

        study_metadata = metadata_file.read_text()
        paper_text = paper_text_file.read_text()
        sample_list = list(set(adata.obs['latch_sample_id']))
        print(f"Harmonizing against sample_list from obs['latch_sample_id']: {sample_list}")

        for var_def in var_defs:
            print(f"Requesting harmonized metadata from model for name: {var_def.name} ; description: {var_def.description}")
            resp = post(
                f"{lcc.nucleus_url}/{lcc.get_harmonized_metadata_endpoint}",
                {
                    "study_metadata": study_metadata,
                    "paper_text": paper_text,
                    "sample_list": sample_list,
                    "var_def": var_to_json(var_def),
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

            annotation_dict[var_def.name] = {"annotations": annotations, "reasoning": reasoning}

    for k, v in annotation_dict.items():
        try:
            adata.obs[k] = adata.obs["latch_sample_id"].map(v["annotations"])
        except Exception as e:
            raise ValueError(f"Issue mapping {k}; {e.with_traceback()}")

    html = build_metadata_report(annotation_dict)
    write_html_report(html, workdir, lcc.harmonize_metadata_report_name)
    if not use_metadata:
        with open(cache_path, "w") as f:
            yaml.safe_dump(annotation_dict, f, default_flow_style=False, indent=2)
        print(f"harmonization metadata written to {cache_path}")
    write_anndata(adata, workdir, lcc.harmonize_metadata_adata_name)

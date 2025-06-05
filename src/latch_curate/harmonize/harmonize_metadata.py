from pathlib import Path
from textwrap import dedent
from html import escape
import json
from typing import Literal, Callable

from anndata import AnnData

from latch_curate.utils import write_html_report, write_anndata
from latch_curate.llm_utils import prompt_model, execute_tool_calls
from latch_curate.harmonize.tools import mondo_search, mondo_get_term, uberon_search, uberon_get_term
from latch_curate.constants import latch_curate_constants as lcc

ControlledMetadataKeys = Literal[
         "latch_subject_id", 
         "latch_condition",
         "latch_disease", 
         "latch_tissue",
         "latch_sample_site",
         "latch_sequencing_platform",
         "latch_organism"
         ]
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

def build_metadata_prompts(study_metadata: str, paper_text: str, sample_list:
                           list[str]) -> dict[ControlledMetadataKeys, str]:
    example = {
        "annotations": {"0": "An example label"},
        "reasoning": "Here is some of thinking:\n",
    }
    output_instruction_snippet = f"""
            Return raw JSON (not markdown) with keys `annotations` and `reasoning`.

            `reasoning` must be a **full markdown document with newlines**, like the example below.

            <example>
            {json.dumps(example)}
            </example>
            """

    # todo(kenny): factor out duplication when these stabilize
    return {
            "latch_subject_id": dedent(f"""
            <sample_list>
            {sample_list}
            </sample_list>

            <paper_text>
            {paper_text}
            </paper_text>

            <study_metadata>
            {study_metadata}
            </study_metadata>

            Create labels to describe unique patient subjects for each of the samples described in <sample_list>.

            Explain your reasoning and show me your chain-of-thought using the information provided in <study_metadata> and <paper_text>

            {output_instruction_snippet}
            """),
            "latch_condition": dedent(f"""
            Now create labels to describe the conditions for each of the samples described in <sample_list>.

            Explain your reasoning and show me your chain-of-thought using the information provided in <study_metadata> and <paper_text>

            Return a dictionary of values with each of <sample_list> values used as the keys.

            {output_instruction_snippet}
            """),
            # todo(kenny): remove
            "latch_disease": dedent(f"""
            <ontology>
            Use the official MONDO ontology ( https://monarchinitiative.org/ ).
            Example of a correct value: "ankylosing spondylitis/MONDO:0005306"
            Use "healthy/" always to refer to healthy.
            </ontology>

            Now create labels to describe the disease for each of the samples in <sample_list>.

            Explain your reasoning using the information provided in <study_metadata> and <paper_text>.

            Using the <ontology>, provide the exact name and ID for each value. Use exactly
            those names and IDs in your final dictionary.

            Show me your sources or chain-of-thought for how you decided on the ID. Then
            confirm you are matching the exact label text from <ontology>. Double check that the
            ontology names and the IDs are correct.

            Return a dictionary of values with each of <sample_list> values used as the
            keys and ontology values as "<name>/<ID>".

            {output_instruction_snippet}
            """),
            "latch_tissue": dedent(f"""
            <ontology>
            Use the official UBERON ontology ( https://obofoundry.org/ontology/uberon.html ).
            Example of a correct value: "blood/UBERON:0000178"
            </ontology>

            Now create labels to describe the tissue for each of the samples in <sample_list>.

            Explain your reasoning using the information provided in <study_metadata> and <paper_text>.

            Using the <ontology>, provide the exact name and ID for each value. Use exactly
            those names and IDs in your final dictionary.

            Show me your sources or chain-of-thought for how you decided on the ID. Then
            confirm you are matching the exact label text from <ontology>. Double check that the
            ontology names and the IDs are correct.

            Return a dictionary of values with each of <sample_list> values used as the
            keys and ontology values as "<name>/<ID>".

            {output_instruction_snippet}
            """),
            "latch_sample_site": dedent(f"""
            <vocabulary>
            lesional, 
            peri-lesional
            normal
            blood
            in vitro
            </vocabulary>

            Now create labels to describe the sampling site for each of the samples in
            <sample_list> using only the exact terms in <vocabulary>.

            Explain your reasoning using the information provided in <study_metadata> and <paper_text>.

            Show me your sources or chain-of-thought for how you decided on the value. then
            confirm you are matching the exact terms in <vocabulary>. double check that the
            terms are correct

            Return a dictionary of values with each of <sample_list> values used as the
            keys.

            {output_instruction_snippet}
            """),
            "latch_sequencing_platform": dedent(f"""
            Now create labels to describe the single cell sequencing platform (usually a
            kit) for each of the samples in <sample_metadata>

            Explain your reasoning using the information provided in <study_metadata> and <paper_text>.

            Show me your sources or chain-of-thought for how you decided on the values.

            {output_instruction_snippet}
            """),
            "latch_organism": dedent(f"""
            Now create a label to describe the species for each sample.

            Explain your reasoning using the information provided in <study_metadata> and <paper_text>.

            Show me your sources or chain-of-thought for how you decided on the values.

            Return a dictionary of values with each of <sample_list> values used as the
            keys.

            {output_instruction_snippet}
            """),
            }

tools_by_key: dict[ControlledMetadataKeys, list[Callable]] = {
         "latch_subject_id": [],
         "latch_condition": [],
         "latch_disease": [mondo_search, mondo_get_term], 
         "latch_tissue": [uberon_search, uberon_get_term], 
         "latch_sample_site": [], 
         "latch_sequencing_platform": [],
         "latch_organism": []
        }

def harmonize_metadata(
    adata: AnnData,
    metadata_file: Path,
    paper_text_file: Path,
    workdir: Path,
):
    workdir.mkdir(exist_ok=True)

    study_metadata = metadata_file.read_text()
    paper_text = paper_text_file.read_text()
    sample_list = set(adata.obs['latch_sample_id'])
    print(f"Harmonizing against sample_list from obs['latch_sample_id']: {sample_list}")

    prompt_dict = build_metadata_prompts(study_metadata, paper_text, sample_list)
    response_dict = {}
    messages = []
    for k, prompt in prompt_dict.items():
        messages.append({"role": "user", "content": prompt_dict[k]})
        while True:
            try:
                print(f"Prompting model for key {k}")
                message_resp_json, tool_calls_resp = prompt_model(messages,
                                                                  tools=tools_by_key[k])

                if tool_calls_resp is not None:
                    tool_out, fn_name, fn_args = execute_tool_calls(tools_by_key[k], tool_calls_resp)
                    messages.append({
                        "role": "user",
                        "content": json.dumps({"tool_out": tool_out, "fn_name": fn_name, "fn_args": fn_args}),
                    })
                    continue

                try:
                    message_resp = json.loads(message_resp_json)
                    annotations: dict = message_resp["annotations"]
                    if set(annotations.keys()) != sample_list:
                        print(f">>> resp: {annotations.keys()}")
                        print("Annotation keys not sample_ids. Trying again.")
                        print(messages)
                        continue

                    reasoning: str = message_resp["reasoning"]
                    print(f"[{k}] annotations >> {annotations}")
                    print(f"[{k}] reasoning >> {reasoning}")
                    response_dict[k] = {"annotations": annotations, "reasoning": reasoning}
                    break
                except Exception:
                    print(f">>> resp: {message_resp_json}")
                    print("Malformed JSON output model. Trying again.")
                    continue

            except Exception as e:
                print(f"dumping error: {e.with_traceback()}")
                print("invalid model response, retrying…")

    for k, v in response_dict.items():
        try:
            adata.obs[k] = adata.obs["latch_sample_id"].map(v["annotations"])
        except Exception as e:
            raise ValueError(f"Issue mapping {k}; {e.with_traceback()}")

    html = build_metadata_report(response_dict)
    write_html_report(html, workdir, lcc.harmonize_metadata_report_name)
    write_anndata(adata, workdir, lcc.harmonize_metadata_adata_name)

from dataclasses import dataclass
from textwrap import dedent
import json
from pathlib import Path
import re

from anndata import AnnData

from latch_curate.llm_utils import prompt_model
from latch_curate.utils import write_html_report
from latch_curate.constants import latch_curate_constants as lcc

@dataclass
class Tag:
    metadata_type: str
    name: str
    obo_id: str


def build_get_paper_info_prompt(paper_text: str) -> str:
    example = {
        "paper_title": "A framework for evaluating ILC2 response in...",
        "paper_abstract": "The role of CD4 T cells in...\n",
    }
    output_instruction_snippet = dedent(
        f"""
        Return raw JSON (not markdown) with keys `paper_title` and `paper_abstract`.

        <example>
        {json.dumps(example)}
        </example>
        """
    )

    return dedent(
        f"""
        <paper_text>
        {paper_text}
        </paper_text>

        Extract the paper title and paper abstract from <paper_text>.

        {output_instruction_snippet}
        """
    )

def build_get_paper_author_contact_info_prompt(paper_text: str) -> str:
    example = {
        "corresponding_author_names": ["John Smith", "Jane Doe"],
        "corresponding_author_emails": ["john_smith@harvard.edu", "jane_doe@harvard.edu"],
    }
    output_instruction_snippet = dedent(
        f"""
        Return raw JSON (not markdown) with keys `corresponding_author_names`
        and `corresponding_author_emails`.

        <example>
        {json.dumps(example)}
        </example>
        """
    )

    return dedent(
        f"""
        <paper_text>
        {paper_text}
        </paper_text>

        Extract the corresponding author names and emails from <paper_text>.

        {output_instruction_snippet}
        """
    )


def build_publish_data(paper_text: str, paper_url: str, gse_id: str, adata: AnnData, workdir: Path, tags: list[Tag]):

    workdir.mkdir(exist_ok=True)

    get_paper_info_prompt = build_get_paper_info_prompt(paper_text)
    get_paper_author_contact_info_prompt = build_get_paper_author_contact_info_prompt(paper_text)

    while True:
        print("Requesting paper metadata from language model")
        message_resp_json, _ = prompt_model([{"role": "user", "content": get_paper_info_prompt}])
        try:
            data = json.loads(message_resp_json)
            paper_title = data["paper_title"]
            paper_abstract = data["paper_abstract"]
            break
        except Exception:
            print("Invalid model response: {message_resp_json}. Trying again")
            continue

    while True:
        print("Requesting corresponding author info from language model")
        message_resp_json, _ = prompt_model([{"role": "user", "content": get_paper_author_contact_info_prompt}])
        try:
            data = json.loads(message_resp_json)
            corresponding_author_names = data["corresponding_author_names"]
            corresponding_author_emails = data["corresponding_author_emails"]
            break
        except Exception:
            print("Invalid model response: {message_resp_json}. Trying again")
            continue

    display_name = f"orion-{gse_id}"

    build_info_file = workdir / lcc.publish_build_info_file_name
    with open(build_info_file, "w") as f:
        data = {
                "info": {
                    "description": paper_abstract,
                    "paper_title": paper_title,
                    "cell_count": adata.n_obs,
                    "paper_url": paper_url,
                    "data_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}",
                    "data_external_id": gse_id,
                    "display_name": display_name,
                    "corresponding_author_names": corresponding_author_names,
                    "corresponding_author_emails": corresponding_author_emails,
                },
        "tags": [{"metadata_type": t.metadata_type, "value": t.name,
                  "ontology_id": t.obo_id} for t in tags],
        }
        json.dump(data, f)
        print(f"Publish build data written to {build_info_file}")

body_pattern = re.compile(r"<body[^>]*>(.*?)</body>", re.S | re.I)

def read_inner_html(path: Path) -> str:
    txt = path.read_text()
    m = body_pattern.search(txt)
    return m.group(1) if m else txt

def build_publish_report(workdir: Path, report_map: dict[str, Path]):
    workdir.mkdir(exist_ok=True)

    tab_buttons = []
    tab_contents = []
    for idx, (name, report_path) in enumerate(report_map.items()):

        assert report_path.exists()

        tab_buttons.append(
            f'<button class="tablinks" onclick="openTab(event, \'tab{idx}\')">{name}</button>'
        )
        tab_contents.append(
            dedent(
                f"""
                <div id="tab{idx}" class="tabcontent" style="overflow:auto;max-height:90vh;">
                    {read_inner_html(report_path)}
                </div>
                """
            )
        )

    html: str = dedent(
        f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="utf-8" />
            <title>latch‑curate – Publish report</title>
            <style>
                body {{ font-family: -apple-system, BlinkMacSystemFont, 'Helvetica', sans-serif; margin:0 }}
                .tabs {{ overflow-x:auto; white-space:nowrap; background:#eee; padding:6px 8px }}
                .tablinks {{ background:#ddd; border:none; outline:none; padding:10px 16px; cursor:pointer; margin-right:4px; border-radius:4px 4px 0 0 }}
                .tablinks:hover {{ background:#ccc }}
                .tablinks.active {{ background:#fff; border-bottom:2px solid #fff }}
                .tabcontent {{ display:none; padding:0 12px 12px 12px; }}
            </style>
            <script>
                function openTab(evt, tabId) {{
                    var i, tabcontent, tablinks;
                    tabcontent = document.getElementsByClassName("tabcontent");
                    for (i = 0; i < tabcontent.length; i++) {{ tabcontent[i].style.display = "none"; }}
                    tablinks = document.getElementsByClassName("tablinks");
                    for (i = 0; i < tablinks.length; i++) {{ tablinks[i].className = tablinks[i].className.replace(" active", ""); }}
                    document.getElementById(tabId).style.display = "block";
                    evt.currentTarget.className += " active";
                }}
                document.addEventListener("DOMContentLoaded", function() {{
                    // auto‑open first available tab
                    var first = document.getElementsByClassName('tablinks')[0];
                    if (first) first.click();
                }});
            </script>
        </head>
        <body>
            <div class="tabs">
                {''.join(tab_buttons)}
            </div>
            {''.join(tab_contents)}
        </body>
        </html>
        """
    )

    write_html_report(html, workdir, lcc.publish_report_name)

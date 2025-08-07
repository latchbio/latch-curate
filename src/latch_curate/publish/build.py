from textwrap import dedent
from pathlib import Path
import re
import yaml

from anndata import AnnData
import requests

from latch_curate.utils import write_html_report
from latch_curate.constants import latch_curate_constants as lcc
from latch_curate.config import user_config
from latch_curate.lint.metadata_validator import validate_harmonized_metadata


def build_publish_data(paper_text: str, paper_url: str, gse_id: str, adata: AnnData, workdir: Path):

    workdir.mkdir(exist_ok=True)

    print("Requesting paper metadata from API")
    resp = requests.post(
        f"{lcc.nucleus_url}/{lcc.get_paper_info_endpoint}",
        json={
            "paper_text": paper_text,
            "session_id": -1
        },
        headers={"Authorization": f"Latch-SDK-Token {user_config.token}"},
        timeout=60
    )
    resp.raise_for_status()
    paper_info_data = resp.json()["data"]
    paper_title = paper_info_data["paper_title"]
    paper_abstract = paper_info_data["paper_abstract"]

    print("Requesting corresponding author info from API")
    resp = requests.post(
        f"{lcc.nucleus_url}/{lcc.get_paper_author_contact_info_endpoint}",
        json={
            "paper_text": paper_text,
            "session_id": -1
        },
        headers={"Authorization": f"Latch-SDK-Token {user_config.token}"},
        timeout=60
    )
    resp.raise_for_status()
    author_info_data = resp.json()["data"]
    corresponding_author_names = author_info_data["corresponding_author_names"]
    corresponding_author_emails = author_info_data["corresponding_author_emails"]

    tags = []
    validation_status = "not_validated"
    schema_path = user_config.metadata_schema_path

    if schema_path.exists():
        print(f"Using metadata schema at: {schema_path}")
        print("Validating harmonized metadata against schema...")

        try:
            validation_result = validate_harmonized_metadata(adata, schema_path)
            tags = validation_result.tags
            validation_status = "passed" if validation_result.is_valid else "failed"

            if not validation_result.is_valid:
                print(f"Metadata validation failed with {len(validation_result.errors)} errors")
        except Exception as e:
            print(f"Error during metadata validation: {e}")
            validation_status = "error"
    else:
        print(f"No metadata schema found at: {schema_path}")
        print("Skipping validation - ensure schema file exists before running publish")

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
                    "corresponding_author_names": corresponding_author_names,
                    "corresponding_author_emails": corresponding_author_emails,
                },
        "validation": {
            "status": validation_status,
            "schema_used": str(schema_path) if schema_path else None,
            "tags_extracted": len(tags),
        },
        "tags": [{"metadata_type": t.metadata_type, "value": t.name,
                  "ontology_id": t.obo_id} for t in tags],
        }
        yaml.safe_dump(data, f, default_flow_style=False, indent=2)
        print(f"Publish build data written to {build_info_file}")
        
    print("\nBuild complete! Please verify the following information:")
    print("=" * 60)
    print(f"Paper Title: {paper_title}")
    print(f"Paper Abstract: {paper_abstract[:200]}...")
    print(f"Cell Count: {adata.n_obs:,}")
    print(f"Authors: {', '.join(corresponding_author_names)}")
    print(f"Email Contacts: {', '.join(corresponding_author_emails)}")
    print(f"Validation Status: {validation_status}")
    print(f"Tags Extracted: {len(tags)}")
    if tags:
        print("Tags:")
        for tag in tags[:5]:
            print(f"  - {tag.metadata_type}: {tag.name}")
        if len(tags) > 5:
            print(f"  ... and {len(tags) - 5} more")
    print("=" * 60)
    print(f"Build file: {build_info_file}")
    print("Next step: latch-curate publish upload")

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

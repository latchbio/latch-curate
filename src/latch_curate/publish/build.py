from dataclasses import dataclass
from textwrap import dedent
import json
from pathlib import Path

from anndata import AnnData

from latch_curate.llm_utils import prompt_model
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


def build_publish_data(paper_text: str, paper_url: str, gse_id: str, adata: AnnData, workdir: Path, tags: list[Tag]):

    workdir.mkdir(exist_ok=True)

    get_paper_info_prompt = build_get_paper_info_prompt(paper_text)

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
                    "display_name": f"orion-{gse_id}"
                },
        "tags": [{"metadata_type": t.metadata_type, "value": t.name,
                  "ontology_id": t.obo_id} for t in tags],
        }
        json.dump(data, f)
        print(f"Publish build data written to {build_info_file}")

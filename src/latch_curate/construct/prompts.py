import json
from textwrap import dedent

def build_construct_counts_instructions(paper_text: str, study_metadata: str):
    return f"""
    <paper_text>
    {paper_text}
    </paper_text>

    <study_metadata>
    {study_metadata}
    </study_metadata>
    """

def build_get_target_cell_count_prompt(paper_text: str, study_metadata: str):

    example = {
        "target_cell_count": 123000,
        "reasoning": "Here is some of thinking:\n",
    }
    output_instruction_snippet = dedent(f"""
    Return raw JSON (not markdown) with keys `target_cell_count` and `reasoning`.

    `reasoning` must be a **full markdown document with newlines**, like the example below.

    <example>
    {json.dumps(example)}
    </example>
    """)

    return dedent(f"""
    {build_construct_counts_instructions(paper_text, study_metadata)}

    Given the <paper_text> and <study_metadata> above, determine the total cell
    count from single cell RNA sequencing.

    {output_instruction_snippet}
    """)


def build_construct_counts_prompt(target_cell_count: int):

    return dedent(f"""
    ## CONTEXT

    You are operating inside a clean working directory that already contains:

    * **Raw data folder**: `data`
    * **Utility library**: `scrna_utils.py`  

    ---

    ## GOAL

    Write a single **driver script** called `build_anndata.py` that:

    1. **Organize** counts and metadata in per-sample folders, extracting or unzipping files as unnecessary.
    2. **Parses** each sample’s count + metadata files using *functions in `scrna_utils.py`* (import and monkey patch as needed)
    3. **Ensures** per-sample AnnData objects satisfy:  
       * `var` index = Ensembl IDs  
       * `var['gene_symbols']` present  
       * all author metadata columns prefixed with `author_`
    4. **Merges** the sample AnnData objects and identify distinct samples with an obs variable named "latch_sample_id". Prefix obs names with sample names to ensure uniqueness.
    5. **Writes** the combined object to `output.h5ad` (do **not** transform or QC the counts).
    6. **Validates** the combined object respects criteria in the ##Validation section.

    ---

    ## GUIDELINES

    * Try to use and monkey patch the helper functions already provided in
    `scrna_utils.py` for relevant components of workflow.
    * Try to match the author cell count provided in <paper_text> with the total cell count in your matrix
    * Try to use all of the relevant information in <paper_text> and <study_metadata> to help you with your task

    ---

    ## VALIDATION


    Validate the following in the constructed AnnData object:

    - there is an obs variable named 'latch_sample_id'
    - the var index is ensembl ids
    - there is a var variable named 'gene_symbols' with the symbols
    - the obs index, var index and var 'gene_symbols' contain unique values'
    - the counts are raw/not transformed (eg. negative, not integers, etc.)
    - the number of rows roughly matches the cell count described in <paper_text>
    - any additional obs variables (eg. `author_` prefixed) have realistic values and not `nan` or similar
    - the count matrix is written to a file named `output.h5ad`
    - there is subject-level metadata available somewhere (if not in `latch_sample_id` then some author variable)
    - the total number of cells is close to {target_cell_count}

    ## STRICT CONSTRAINTS

    *  **Do not** normalise, log-transform, or filter the counts.  

    After you finish writing build_anndata.py, execute it with "/Users/kenny/latch/latch-curate/construct/.venv/bin/python3 build_anndata.py" and do not exit until the file output.h5ad exists.
    """)

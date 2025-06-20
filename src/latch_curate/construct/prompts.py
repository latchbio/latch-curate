import re
from textwrap import dedent


def build_review_prompt(paper_text: str, study_metadata: str, driver_script:
                        str, agent_prompt: str, agent_logs: str, query: str):
    return f"""
    <paper_text>
    {paper_text}
    </paper_text>

    <study_metadata>
    {study_metadata}
    </study_metadata>

    <driver_script>
    {driver_script}
    </driver_script>

    <agent_prompt>
    {agent_prompt}
    </agent_prompt>

    <agent_logs>
    {agent_logs}
    </agent_logs>

    Use all of the context above to answer the following question about past agent
    behavior:

    {query}
    """

def build_construct_counts_instructions(paper_text: str, study_metadata: str):
    return f"""
    <paper_text>
    {paper_text}
    </paper_text>

    <study_metadata>
    {study_metadata}
    </study_metadata>
    """

validation_failure_pattern = r"<validation_failure>.*?</validation_failure>"

def add_or_replace_validation_failure(prompt: str, validation_failure: str):
    tagged_failure = dedent(f"""\
    <validation_failure>
    {validation_failure.rstrip()}
    </validation_failure>
    """)
    if re.search(validation_failure_pattern, prompt):
         return re.sub(
             validation_failure_pattern,
             tagged_failure,
             prompt,
         )
    return f"{prompt}\n\n{tagged_failure}"

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
    * The paper reports a total cell count around {target_cell_count}. Use this
    as a guideline, especially if you are below this value, as it indicates you
    should look for more data. NEVER downsample or drop raw data if you are
    above this target value.
    * Prioritize the information in <validation_failure>, if it exists. It
    contains failed validation tests that check the criteria in ## VALIDATION
    from your last attempt to perform this task.
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

    ## STRICT CONSTRAINTS

    * **Do not** normalise, log-transform, or filter the counts.  
    * Under **no circumstances** should the script subset, sample, or
    downsample the raw data.  It must load every barcode and every nonzero
    count value.

    After you finish writing build_anndata.py, execute it with "/Users/kenny/latch/latch-curate/construct/.venv/bin/python3 build_anndata.py" and do not exit until the file output.h5ad exists.
    """)

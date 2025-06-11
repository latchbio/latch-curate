from urllib.parse import quote
from latch_curate.tinyrequests import get

OLS_BASE = "https://www.ebi.ac.uk/ols4/api"


def _fetch_json(url: str) -> dict:
    resp = get(url)
    resp.raise_for_status()
    return resp.json()


def _ols_search(term: str, *, ontology: str, rows: int = 20) -> dict:
    q = quote(term, safe="")
    url = f"{OLS_BASE}/search?q={q}&ontology={ontology}&rows={rows}"
    return _fetch_json(url)


def _ols_get_term(curie: str, *, ontology: str) -> dict:
    curie_us = curie.replace(":", "_")
    iri_enc = quote(f"http://purl.obolibrary.org/obo/{curie_us}", safe="")
    url = f"{OLS_BASE}/ontologies/{ontology}/terms/{iri_enc}"
    return _fetch_json(url)


def mondo_search(term: str, rows: int = 20) -> dict:
    """
    Search the MONDO disease ontology with a free-text query.

    Intended for use by an agent tool that needs to map arbitrary disease
    names or synonyms to MONDO concepts.  The helper issues a GET request to
    the EBI OLS `/search` endpoint restricted to the **MONDO** ontology and
    returns the raw JSON payload.

    Parameters
    ----------
    term : str
        Free-text query (e.g. ``"breast cancer"``, ``"glioblastoma"``,
        ``"AML"``).
    rows : int, default 20
        Maximum number of search hits to return.

    Returns
    -------
    dict
        Parsed OLS search result.  Relevant hits are in
        ``result["response"]["docs"]``.

    Example
    -------
    >>> hit = mondo_search("glioblastoma")["response"]["docs"][0]
    >>> hit["curie"], hit["label"]
    ('MONDO:0007901', 'glioblastoma')
    """
    return _ols_search(term, ontology="mondo", rows=rows)


def mondo_get_term(mondo_id: str) -> dict:
    """
    Retrieve a single MONDO concept by its CURIE.

    Designed for agent tool workflows that already resolved a MONDO identifier
    (e.g. via :pyfunc:`mondo_search`) and now need full term metadata—
    definition, synonyms, parents, children, xrefs, etc.

    Parameters
    ----------
    mondo_id : str
        MONDO compact identifier such as ``"MONDO:0007901"``.

    Returns
    -------
    dict
        Complete OLS JSON record for the requested MONDO term.

    Example
    -------
    >>> mondo_get_term("MONDO:0007901")["annotation"]["definition"][0]
    'A malignant, aggressive astrocytoma…'
    """
    return _ols_get_term(mondo_id, ontology="mondo")


def uberon_search(term: str, rows: int = 20) -> dict:
    """
    Search the UBERON multi-species anatomy ontology via OLS.

    Useful for agents that must map anatomical phrases to UBERON IDs before
    downstream processing (e.g. tissue harmonisation in single-cell metadata).

    Parameters
    ----------
    term : str
        Free-text anatomy query (e.g. ``"liver"``, ``"left kidney"``,
        ``"cardiac muscle"``).
    rows : int, default 20
        Maximum number of hits to return.

    Returns
    -------
    dict
        Raw OLS search response with matching UBERON terms ranked by score.

    Example
    -------
    >>> uberon_search("liver")["response"]["docs"][0]["curie"]
    'UBERON:0002107'
    """
    return _ols_search(term, ontology="uberon", rows=rows)


def uberon_get_term(uberon_id: str) -> dict:
    """
    Fetch detailed metadata for a single UBERON term.

    Parameters
    ----------
    uberon_id : str
        UBERON compact identifier (e.g. ``"UBERON:0002107"``).

    Returns
    -------
    dict
        Full JSON record for the anatomy concept, including label, definition,
        synonyms, and hierarchical relationships.

    Example
    -------
    >>> uberon_get_term("UBERON:0002107")["label"]
    'liver'
    """
    return _ols_get_term(uberon_id, ontology="uberon")

def efo_search(term: str, rows: int = 20) -> dict:
    """Search the **EFO** (Experimental Factor Ontology) by free‑text query.

    EFO covers human traits, diseases, sample characteristics, and experimental
    variables. These helpers are useful for mapping arbitrary user‑supplied
    labels to stable EFO identifiers prior to downstream harmonisation.

    Parameters
    ----------
    term : str
        Free‑text query (e.g. ``"cancer"``, ``"type 2 diabetes"``).
    rows : int, default 20
        Maximum number of hits to return.

    Returns
    -------
    dict
        Parsed OLS search result. Hits are found under
        ``result["response"]["docs"]``.

    Example
    -------
    >>> hit = efo_search("breast carcinoma")["response"]["docs"][0]
    >>> hit["curie"], hit["label"]
    ('EFO:0000305', 'breast carcinoma')
    """
    return _ols_search(term, ontology="efo", rows=rows)


def efo_get_term(efo_id: str) -> dict:
    """Fetch detailed metadata for a single EFO term by CURIE.

    Parameters
    ----------
    efo_id : str
        Compact identifier such as ``"EFO:0000305"``.

    Returns
    -------
    dict
        Complete OLS JSON record, including label, definition, synonyms, and
        hierarchical relationships.

    Example
    -------
    >>> efo_get_term("EFO:0000305")["label"]
    'breast carcinoma'
    """
    return _ols_get_term(efo_id, ontology="efo")

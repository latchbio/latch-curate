import urllib.parse
import requests

from . import common

EFO_BASE_URL = "https://www.ebi.ac.uk/ols/api/ontologies/efo/terms"

def efo_id_to_iri(efo_id: str) -> str:
    prefix, local = efo_id.split(":", 1)
    return f"http://www.ebi.ac.uk/efo/{prefix}_{local}"

def validate_efo(val: str) -> bool:
    if val == common.unknown_val:
        return True

    if "/" not in val:
        return False

    name_part, efo_id = val.split("/", 1)
    name_part = name_part.strip()

    if not efo_id.startswith("EFO:"):
        return False

    iri = efo_id_to_iri(efo_id)
    query_url = f"{EFO_BASE_URL}?iri={urllib.parse.quote(iri, safe='')}"

    try:
        resp = requests.get(query_url, timeout=10)
        resp.raise_for_status()
        data = resp.json()
    except requests.RequestException:
        return False

    page = data.get("page", {})
    if page.get("totalElements", 0) == 0:
        return False

    terms = data.get("_embedded", {}).get("terms", [])
    if not terms:
        return False

    term = terms[0]
    label = term.get("label", "") or ""

    return name_part.strip().lower() == label.lower()

import urllib.parse

import requests

from . import common

OLS_BASE_URL = "https://www.ebi.ac.uk/ols/api/ontologies/mondo/terms"

healthy_val = "healthy/"


def mondo_id_to_iri(mondo_id: str) -> str:
    underscore_id = mondo_id.replace(":", "_", 1)
    return f"http://purl.obolibrary.org/obo/{underscore_id}"


def validate_mondo(val: str) -> bool:
    if val == common.unknown_val:
        return True
    if val == healthy_val:
        return True
    if "/" not in val:
        return False
    name_part, mondo_id = val.split("/", 1)
    name_part = name_part.strip()

    if not mondo_id.startswith("MONDO:"):
        return False

    iri = mondo_id_to_iri(mondo_id)
    query_url = f"{OLS_BASE_URL}?iri={urllib.parse.quote(iri, safe='')}"

    try:
        resp = requests.get(query_url, timeout=10)
        if resp.status_code != 200:
            return False

        data = resp.json()
    except requests.RequestException:
        return False

    if "page" not in data or data["page"].get("totalElements", 0) == 0:
        return False

    terms_list = data["_embedded"]["terms"]
    if not terms_list:
        return False

    term = terms_list[0]
    label = term.get("label", "") or ""
    synonyms = term.get("synonyms", [])

    name_lc = name_part.lower()
    label_lc = label.lower()

    if name_lc == label_lc:
        return True

    return False

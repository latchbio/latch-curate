import urllib.parse

import requests

UBERON_BASE_URL = "https://www.ebi.ac.uk/ols/api/ontologies/uberon/terms"

def uberon_id_to_iri(uberon_id: str) -> str:
    underscore_id = uberon_id.replace(":", "_", 1)
    return f"http://purl.obolibrary.org/obo/{underscore_id}"

def validate_uberon(val: str) -> bool:

    if "/" not in val:
        return False
    label_part, uberon_id = val.split("/", 1)
    label_part = label_part.strip()
    if not uberon_id.startswith("UBERON:"):
        return False

    iri = uberon_id_to_iri(uberon_id)

    query_url = f"{UBERON_BASE_URL}?iri={urllib.parse.quote(iri, safe='')}"

    try:
        resp = requests.get(query_url, timeout=10)
        if resp.status_code != 200:
            return False
        data = resp.json()
    except requests.RequestException:
        return False

    page_info = data.get("page", {})
    if page_info.get("totalElements", 0) == 0:
        return False

    terms_list = data["_embedded"].get("terms", [])
    if not terms_list:
        return False

    term = terms_list[0]
    label = term.get("label", "") or ""

    label_lc = label.lower()
    user_label_lc = label_part.lower()

    if user_label_lc == label_lc:
        return True

    return False

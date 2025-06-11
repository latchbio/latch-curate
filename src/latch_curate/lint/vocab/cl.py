import urllib.parse

import requests

from . import common

CL_BASE_URL = "https://www.ebi.ac.uk/ols/api/ontologies/cl/terms"


def cl_id_to_iri(cl_id: str) -> str:
    underscore_id = cl_id.replace(":", "_", 1)
    return f"http://purl.obolibrary.org/obo/{underscore_id}"


def validate_cl(user_val: str) -> bool:
    if user_val == common.unknown_val:
        return True

    if "/" not in user_val:
        return False
    label_part, cl_id = user_val.split("/", 1)
    label_part = label_part.strip()
    if not cl_id.startswith("CL:"):
        return False

    iri = cl_id_to_iri(cl_id)

    query_url = f"{CL_BASE_URL}?iri={urllib.parse.quote(iri, safe='')}"
    try:
        resp = requests.get(query_url, timeout=10)
        if resp.status_code != 200:
            return False  # If 404, 400, etc., it's not recognized
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
    official_label = term.get("label", "") or ""

    user_label_lc = label_part.lower()
    official_label_lc = official_label.lower()

    if user_label_lc == official_label_lc:
        return True

    return False

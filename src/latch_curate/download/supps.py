import requests
from lxml import html
from pathlib import Path

def download_gse_supps(gse_id: str, target_dir: Path):

    series_dir = f"{gse_id[:-3]}nnn"
    root = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{series_dir}/{gse_id}/suppl/"

    resp = requests.get(root)
    resp.raise_for_status()
    tree = html.fromstring(resp.content)
    hrefs = tree.xpath("//a/@href")

    # todo(kenny): kind of bad
    # drop parent‐dir link and any policy redirects
    files = [h for h in hrefs if not h.endswith("/") and "vulnerability" not in h]

    target_dir.mkdir(exist_ok=True)
    for fname in files:
        url = root + fname
        dest = target_dir / fname

        print(f"> Downloading {fname} ...")
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(dest, "wb") as f:
                for chunk in r.iter_content(chunk_size=1024*1024):
                    if chunk:
                        f.write(chunk)

    print(f"All files for {gse_id} downloaded into {target_dir!r}.")

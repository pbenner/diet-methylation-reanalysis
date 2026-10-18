"""Fetch the public processed dataset; never download the 192 raw IDAT files."""
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
ACCESSION = "E-MTAB-17025"
API = f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{ACCESSION}"


def download(url, path, expected_size=None):
    if path.exists() and (expected_size is None or path.stat().st_size == expected_size):
        return
    temporary = path.with_suffix(path.suffix + ".part")
    print(f"Downloading {url}", flush=True)
    with urlopen(url, timeout=120) as response, temporary.open("wb") as output:
        while chunk := response.read(8*1024*1024):
            output.write(chunk)
    if expected_size is not None and temporary.stat().st_size != expected_size:
        raise ValueError(f"Unexpected file size: {temporary}")
    temporary.replace(path)


def files_in(node):
    if isinstance(node, dict):
        if "path" in node and "size" in node:
            yield node
        for value in node.values():
            yield from files_in(value)
    elif isinstance(node, list):
        for value in node:
            yield from files_in(value)


def main():
    raw = ROOT / "data/raw"
    raw.mkdir(parents=True, exist_ok=True)
    study_path, info_path = raw/f"{ACCESSION}.json", raw/f"{ACCESSION}-info.json"
    download(API, study_path)
    download(API+"/info", info_path)
    study, info = json.loads(study_path.read_text()), json.loads(info_path.read_text())
    entries = {f["path"]: f for f in files_in(study)}
    manifest = []
    for name in (f"{ACCESSION}.sdrf.txt", f"{ACCESSION}.idf.txt", "b_values4.txt"):
        url = info["httpLink"] + "/Files/" + name
        path = raw/name
        download(url, path, entries[name]["size"])
        with path.open("rb") as handle:
            checksum = hashlib.file_digest(handle, "sha256").hexdigest()
        manifest.append({"file": name, "url": url, "bytes": path.stat().st_size, "sha256": checksum})
    (raw/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")


if __name__ == "__main__":
    main()

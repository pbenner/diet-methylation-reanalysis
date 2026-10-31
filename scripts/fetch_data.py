"""Fetch the public processed dataset; never download the 192 raw IDAT files."""
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
ACCESSION = "E-MTAB-17025"
API = f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{ACCESSION}"


def checksum(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def verify_hash(path, expected_hash):
    if expected_hash is not None and checksum(path) != expected_hash:
        raise ValueError(f'{path}: SHA-256 differs from the archived input. '
                         'Keep annotation versions separate; do not silently replace the pinned data.')


def download(url, path, expected_size=None, expected_hash=None):
    if path.exists() and (expected_size is None or path.stat().st_size == expected_size):
        verify_hash(path, expected_hash)
        return
    temporary = path.with_suffix(path.suffix + ".part")
    print(f"Downloading {url}", flush=True)
    with urlopen(url, timeout=120) as response, temporary.open("wb") as output:
        while chunk := response.read(8*1024*1024):
            output.write(chunk)
    if expected_size is not None and temporary.stat().st_size != expected_size:
        raise ValueError(f"Unexpected file size: {temporary}")
    verify_hash(temporary, expected_hash)
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
    study = json.loads(study_path.read_text())
    entries = {f["path"]: f for f in files_in(study)}
    manifest_path = raw/'manifest.json'
    archived = {e['file']:e for e in json.loads(manifest_path.read_text())} if manifest_path.exists() else {}
    manifest = []
    for name in (f"{ACCESSION}.sdrf.txt", f"{ACCESSION}.idf.txt", "b_values4.txt"):
        # The stable download endpoint follows BioStudies storage migrations.
        url = f'https://www.ebi.ac.uk/biostudies/files/{ACCESSION}/{name}'
        path = raw/name
        pinned = archived.get(name, {})
        download(url, path, pinned.get('bytes',entries[name]['size']), pinned.get('sha256'))
        manifest.append({"file": name, "url": url, "bytes": path.stat().st_size, "sha256": checksum(path)})
    # Preserve the original provenance record, including its historical URLs.
    if not archived:
        manifest_path.write_text(json.dumps(manifest, indent=2)+"\n")


if __name__ == "__main__":
    main()

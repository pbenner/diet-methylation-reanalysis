"""Validate deposited annotations and convert the unusual /t-delimited beta matrix."""
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def read_samples(path):
    samples = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample = {
                "sample_id": row["Source Name"],
                "array_id": "ID_" + row["Array Data File"].rsplit("_", 1)[0],
                "participant": row["Characteristics[individual]"],
                "time": int(row["Factor Value[time]"]),
                "diet_code": row["Characteristics[diet]"],
                "diet_named": row["Characteristics[disease]"],
                "age": float(row["Characteristics[age]"]),
                "sex": row["Characteristics[sex]"],
            }
            if sample["diet_named"] != row["Factor Value[disease]"]:
                raise ValueError("Conflicting named diet fields")
            key = sample["sample_id"]
            if key in samples and samples[key] != sample:
                raise ValueError(f"Red/green metadata disagree: {key}")
            samples[key] = sample
    if len({s["array_id"] for s in samples.values()}) != len(samples):
        raise ValueError("Duplicate array IDs")
    participants = defaultdict(list)
    for sample in samples.values():
        participants[sample["participant"]].append(sample)
    for person, pair in participants.items():
        if len(pair) != 2 or {s["time"] for s in pair} != {0, 1}:
            raise ValueError(f"Incomplete/duplicated pair: {person}")
        for field in ("diet_code", "diet_named", "sex", "age"):
            if len({s[field] for s in pair}) != 1:
                raise ValueError(f"Inconsistent {field}: {person}")
    return list(samples.values())


def main():
    raw, out = ROOT / "data/raw", ROOT / "data/processed"
    out.mkdir(parents=True, exist_ok=True)
    samples = read_samples(raw / "E-MTAB-17025.sdrf.txt")
    source = raw / "b_values4.txt"
    expected_size = 1494040132
    if source.stat().st_size != expected_size:
        raise ValueError("Matrix size does not match deposited file; download may be incomplete")
    digest = hashlib.sha256()
    seen = set()
    lower, upper = np.inf, -np.inf
    with source.open("rb") as handle, (out / "beta.bin.part").open("wb") as binary, (out / "cpg_ids.txt.part").open("w") as ids:
        header = handle.readline()
        digest.update(header)
        array_ids = [x.strip('"') for x in header.decode().strip().split("/t")]
        by_array = {s["array_id"]: s for s in samples}
        if len(set(array_ids)) != len(array_ids) or set(array_ids) != set(by_array):
            raise ValueError("Matrix and annotation sample IDs do not match exactly")
        samples = [by_array[key] for key in array_ids]
        for count, line in enumerate(handle, 1):
            digest.update(line)
            cpg, values = line.decode().strip().split("/t", 1)
            cpg = cpg.strip('"')
            if cpg in seen:
                raise ValueError(f"Duplicated CpG: {cpg}")
            seen.add(cpg)
            values = np.fromstring(values.replace("/t", " "), sep=" ", dtype="<f8")
            if values.size != len(samples) or not np.isfinite(values).all():
                raise ValueError(f"Malformed/nonfinite row: {cpg}")
            if not ((values > 0) & (values < 1)).all():
                raise ValueError(f"Beta values outside the open interval (0, 1): {cpg}")
            lower, upper = min(lower, float(values.min())), max(upper, float(values.max()))
            binary.write(values.tobytes())
            ids.write(cpg + "\n")
            if count % 100000 == 0:
                print(f"Validated {count:,} CpGs", flush=True)
    if count != 812934:
        raise ValueError(f"Expected 812934 CpGs from paper; found {count}")
    (out / "beta.bin.part").replace(out / "beta.bin")
    (out / "cpg_ids.txt.part").replace(out / "cpg_ids.txt")
    with (out / "samples.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(samples[0]))
        writer.writeheader()
        writer.writerows(samples)
    baseline = [s for s in samples if s["time"] == 0]
    audit = {
        "matrix_shape": [count, len(samples)], "dtype": "little-endian float64, C row order",
        "beta_min": lower, "beta_max": upper, "matrix_sha256": digest.hexdigest(),
        "participants": len(baseline),
        "numeric_group_counts": dict(Counter(s["diet_code"] for s in baseline)),
        "named_group_counts": dict(Counter(s["diet_named"] for s in baseline)),
        "cross_tab": {f"{code}|{name}": n for (code, name), n in Counter((s["diet_code"], s["diet_named"]) for s in baseline).items()},
        "status": "Conflicting diet annotations. Numeric codes have no verified diet-name mapping; all analyses provisional.",
    }
    (out / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()

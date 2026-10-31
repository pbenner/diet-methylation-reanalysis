"""Apply four explicit visit corrections without changing the deposited inputs."""
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def apply_visit_corrections(samples, corrections):
    result = samples.copy()
    if samples.array_id.duplicated().any() or corrections.array_id.duplicated().any():
        raise ValueError('Array IDs must be unique')
    for row in corrections.itertuples(index=False):
        selected = samples.array_id.eq(row.array_id)
        if selected.sum() != 1:
            raise ValueError(f'Expected exactly one array: {row.array_id}')
        sample = samples.loc[selected].iloc[0]
        if str(sample.participant) != str(row.participant) or sample.time != row.original_time:
            raise ValueError(f'Original participant/time differs for {row.array_id}; do not apply twice')
        if row.original_time not in (0, 1) or row.corrected_time != 1-row.original_time:
            raise ValueError('Expected a reversal of visits 0/1')
        result.loc[selected, 'time'] = row.corrected_time
    if not result.groupby('participant').time.apply(lambda x: sorted(x)==[0, 1]).all():
        raise ValueError('Corrections must preserve one sample per visit per participant')
    pd.testing.assert_frame_equal(samples.drop(columns='time'), result.drop(columns='time'))
    return result


def main():
    source = ROOT/'data/processed/samples.csv'
    patch = ROOT/'data/corrections/visit_labels.csv'
    destination = ROOT/'data/processed/samples_corrected_visits.csv'
    samples = pd.read_csv(source, dtype={'participant':str, 'diet_code':str})
    corrections = pd.read_csv(patch, dtype={'participant':str})
    if len(corrections)!=4 or set(corrections.participant)!={'5','8'}:
        raise ValueError('Expected exactly the four documented corrections for donors 5 and 8')
    result = apply_visit_corrections(samples, corrections)
    result.to_csv(destination, index=False)
    report = {
        'scenario':'corrected_visits', 'corrected_rows':int((samples.time!=result.time).sum()),
        'only_time_changed':True, 'sample_and_matrix_column_order_preserved':True,
        'source_samples_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'correction_sha256':hashlib.sha256(patch.read_bytes()).hexdigest(),
        'corrected_samples_sha256':hashlib.sha256(destination.read_bytes()).hexdigest(),
        'corrections':corrections.to_dict(orient='records'),
        'provenance':'data/corrections/provenance.json',
    }
    (ROOT/'data/processed/visit_correction.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()

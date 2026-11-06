"""Extract the author-supplied workbook into a minimal, validated allocation table."""
import argparse
import hashlib
import json
import re
from math import comb, prod
from pathlib import Path
from zipfile import ZipFile
from xml.etree import ElementTree as ET

import pandas as pd

from randomization import align_blocks

ROOT = Path(__file__).resolve().parents[1]
NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}


def read_workbook(path):
    # Read stored cell values using the standard library; no Excel macros or formulas run.
    with ZipFile(path) as archive:
        strings = [''.join(item.itertext()) for item in
                   ET.fromstring(archive.read('xl/sharedStrings.xml')).findall('m:si', NS)]
        sheet = ET.fromstring(archive.read('xl/worksheets/sheet1.xml'))
    rows = []
    for row in sheet.findall('.//m:row', NS):
        values = {}
        for cell in row.findall('m:c', NS):
            if cell.find('m:f', NS) is not None:
                raise ValueError('Unexpected formula in allocation worksheet')
            value = cell.find('m:v', NS)
            if value is not None:
                values[re.sub(r'\d+$', '', cell.attrib['r'])] = (
                    strings[int(value.text)] if cell.attrib.get('t') == 's' else value.text)
        if row.attrib['r'] == '1':
            if values != {'A': 'PROB-ID', 'C': 'diet', 'D': 'block'}:
                raise ValueError('Unexpected worksheet headings')
            continue
        if not values:
            continue
        if set(values) != {'A', 'B', 'C', 'D'}:
            raise ValueError('Incomplete allocation row')
        person = str(int(values['B']))
        if int(values['A']) != 142000 + int(person):
            raise ValueError('Full study ID and short participant ID disagree')
        match = re.fullmatch(r'([1-4]) - (normal|replacement block)', values['D'])
        if not match or values['C'] not in ('vegan', 'meat-rich'):
            raise ValueError('Unrecognized block or diet')
        number, kind = match.groups()
        kind = 'normal' if kind == 'normal' else 'replacement'
        rows.append({'participant': person, 'diet_code': '0' if values['C'] == 'vegan' else '1',
                     'diet': values['C'], 'block': f'{kind}_{number}',
                     'planned_block_size': 13 if kind == 'normal' else 4})
    return pd.DataFrame(rows).sort_values('participant', key=lambda col: col.astype(int)).reset_index(drop=True)


def main(source):
    table = read_workbook(source)
    samples = pd.read_csv(ROOT/'data/processed/samples_corrected_visits.csv',
                          dtype={'participant': str, 'diet_code': str})
    baseline = samples[samples.time == 0].set_index('participant').sort_index()
    align_blocks(table, baseline.index, baseline.diet_code.eq('0').to_numpy())
    if len(table) != 48 or table.diet_code.value_counts().to_dict() != {'0': 24, '1': 24}:
        raise ValueError('Unexpected participant or diet totals')
    counts = pd.crosstab(table.block, table.diet_code).reindex(columns=['0', '1'], fill_value=0)
    expected_blocks = {f'normal_{i}' for i in range(1, 5)} | {'replacement_1', 'replacement_2'}
    if set(counts.index) != expected_blocks:
        raise ValueError('Expected four normal and two replacement blocks')
    out = ROOT/'data/randomization'
    out.mkdir(exist_ok=True)
    destination = out/'participant_blocks.csv'
    table.to_csv(destination, index=False)
    report = {
        'recorded_on': '2026-11-06', 'source': source.name,
        'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'source_description': 'Allocation workbook supplied by Max Storz to the project owner.',
        'worksheet': 'Tabelle1',
        'participant_mapping': 'Column B matches all 48 deposited participant IDs; column A equals 142000 plus column B.',
        'author_clarification': '0 = vegan, 1 = meat-rich; external statistician created the allocation list; planned blocks of 13 through observation 52, then replacement blocks of 4.',
        'output_file': str(destination.relative_to(ROOT)),
        'output_sha256': hashlib.sha256(destination.read_bytes()).hexdigest(),
        'all_48_numeric_diet_assignments_match': True,
        'planned_block_sizes': {'normal': 13, 'replacement': 4},
        'analyzed_block_counts': {b: {'vegan': int(r['0']), 'meat_rich': int(r['1'])} for b, r in counts.iterrows()},
        'conditional_assignment_space': prod(comb(int(r.sum()), int(r['0'])) for _, r in counts.iterrows()),
        'conditioning': 'Shuffle labels within the supplied blocks, preserving observed counts among analyzed participants in each block; singleton replacement_2 remains fixed.',
        'remaining_limit': 'The workbook covers the analyzed subset. Conditional validity requires exchangeability after subset selection; the full recruitment, replacement and missingness mechanism is not reconstructed.',
    }
    (out/'provenance.json').write_text(json.dumps(report, indent=2)+'\n')
    print(counts.rename(columns={'0': 'vegan', '1': 'meat-rich'}).to_string())
    print(f'Validated all 48 participants; saved {destination}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=ROOT/'ID_BlockAssignment_DrBenner.xlsx')
    main(parser.parse_args().source)

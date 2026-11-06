"""Participant-level allocation restrictions, shared by simulation and examples."""
import hashlib

import numpy as np
import pandas as pd

BLOCK_SCENARIO = 'blocked_randomization'
BLOCK_FILE = 'data/randomization/participant_blocks.csv'


def align_blocks(table, people, labels):
    """Match explicit participant IDs and verify every observed diet assignment."""
    people = list(people)
    labels = np.asarray(labels)
    if labels.dtype != bool or labels.shape != (len(people),):
        raise ValueError('Expected one boolean diet label per participant')
    if (len(set(people)) != len(people) or table.participant.isna().any()
            or table.participant.duplicated().any()
            or set(table.participant) != set(people)):
        raise ValueError('Block table must identify exactly the analyzed participants once')
    aligned = table.set_index('participant').loc[people]
    if (aligned.block.isna().any() or aligned.block.eq('').any()
            or not aligned.diet_code.isin(['0', '1']).all()):
        raise ValueError('Missing block or invalid diet code')
    if not np.array_equal(aligned.diet_code.eq('0').to_numpy(), labels):
        raise ValueError('Author diet assignments disagree with analyzed numeric codes')
    blocks = aligned.block.to_numpy()
    return tuple(np.flatnonzero(blocks == b) for b in sorted(set(blocks)))


def allocation(root, people, labels, scenario, field='diet_code'):
    if scenario != BLOCK_SCENARIO:
        return None, {'permutation_scheme': 'unrestricted',
                      'assumption': 'unrestricted participant exchangeability conditional on total group counts'}
    if field != 'diet_code':
        raise ValueError('Block restrictions require the author-verified numeric diet codes')
    path = root / BLOCK_FILE
    table = pd.read_csv(path, dtype=str)
    indices = align_blocks(table, people, labels)
    counts = pd.crosstab(table.block, table.diet_code).reindex(columns=['0', '1'], fill_value=0)
    return indices, {
        'permutation_scheme': 'within_randomization_blocks',
        'block_file': BLOCK_FILE,
        'block_file_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'block_counts': {block: {'vegan': int(row['0']), 'meat_rich': int(row['1'])}
                         for block, row in counts.iterrows()},
        'diet_code_mapping': {'0': 'vegan', '1': 'meat-rich'},
        'assumption': ('joint exchangeability within the supplied randomization blocks, conditional '
                       'on the analyzed participants and observed diet counts in each block; '
                       'requires the analyzed-subset selection to permit this exchangeability'),
    }


def draw_assignment(rng, labels, block_indices=None):
    """Draw uniformly within each block; retain the historical RNG stream if unrestricted."""
    if block_indices is None:
        return rng.permutation(labels)
    assignment = labels.copy()
    for indices in block_indices:
        assignment[indices] = rng.permutation(labels[indices])
    return assignment

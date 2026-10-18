"""Check source alignment and saved results, with independent real-data R fits.

Requires the downloaded matrix, prepared data, full limma tables and local R.
Writes audit evidence only; never replaces the original analysis results.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def validate_saved_summaries(folder):
    """Cross-check every saved count summary against limma counts and paired rows."""
    def require(condition, message):
        if not condition:
            raise ValueError(message)

    def compare(row, expected, label):
        for key, value in expected.items():
            require(key in row and np.isfinite(row[key])
                    and np.isclose(row[key], value, rtol=1e-12, atol=1e-12),
                    f'{label}: inconsistent {key}')

    observed = pd.read_csv(folder/'limma_summary.csv').set_index('time')
    metadata = json.loads((folder/'permutation_run.json').read_text())
    counts = pd.read_csv(folder/'permutation_null_counts.csv')
    require(observed.index.is_unique and set(observed.index)=={'before','after'},
            'Expected one observed count at each visit')
    require(list(counts.columns)==['before','after'] and len(counts)>0,
            'Expected paired before/after count columns')
    values = counts.to_numpy(dtype=float)
    require(np.isfinite(values).all() and (values>=0).all()
            and (values==np.floor(values)).all(), 'Invalid shuffled counts')
    require(metadata['group_field']=='diet_code' and len(counts)==metadata['repetitions'],
            'Count distribution does not match run metadata')
    for time in ('before','after'):
        value = observed.loc[time,'highlighted']
        require(np.isfinite(value) and value>=0 and value==int(value), 'Invalid observed count')
        require(value==metadata['observed_counts'][time], 'Observed counts differ from run metadata')

    changes = counts.after-counts.before
    paired = pd.read_csv(folder/'permutation_null_changes.csv')
    require(np.array_equal(paired['shuffle'],np.arange(1,len(counts)+1)),
            'Shuffle identifiers are not consecutive unique assignments')
    require(np.array_equal(paired[['before','after']],counts)
            and np.array_equal(paired.change,changes), 'Saved change rows are not aligned')
    change_table = pd.read_csv(folder/'change_test_summary.csv')
    require(len(change_table)==1, 'Expected one change-test summary')
    d = int(observed.loc['after','highlighted']-observed.loc['before','highlighted'])
    extreme = int((changes>=d).sum())
    p = (extreme+1)/(len(counts)+1)
    compare(change_table.iloc[0], {
        'observed_before':observed.loc['before','highlighted'],
        'observed_after':observed.loc['after','highlighted'], 'observed_change':d,
        'permutations':len(counts), 'shuffles_ge_observed':extreme,
        'shuffled_median_change':changes.median(),
        'shuffled_q025_change':changes.quantile(.025),
        'shuffled_q975_change':changes.quantile(.975),
        'change_p_upper':p, 'change_p_mc_se_approx':np.sqrt(p*(1-p)/len(counts)),
    }, 'Change summary')
    visit_summary = pd.read_csv(folder/'count_test_summary.csv').set_index('time')
    require(visit_summary.index.is_unique and set(visit_summary.index)=={'before','after'},
            'Expected one count-test summary per visit')
    for time in ('before','after'):
        extreme_v = int((counts[time]>=observed.loc[time,'highlighted']).sum())
        pv = (extreme_v+1)/(len(counts)+1)
        compare(visit_summary.loc[time], {
            'observed_count':observed.loc[time,'highlighted'], 'permutations':len(counts),
            'shuffles_ge_observed':extreme_v, 'count_p_upper':pv,
            'shuffled_median':counts[time].median(),
            'shuffled_q025':counts[time].quantile(.025),
            'shuffled_q975':counts[time].quantile(.975),
            'count_p_mc_se_approx':np.sqrt(pv*(1-pv)/len(counts)),
        }, time+' count summary')
    return {'all_summary_fields_and_shuffle_ids_match':True,
            'observed_change':d,'permutations':len(counts),
            'shuffles_ge_observed':extreme,'p_upper':p}


def replay_permutations(seed=None):
    """Replay every saved shuffle using the current code, writing only audit output."""
    import permutation

    started = time.monotonic()
    original = ROOT/'results/diet_code/permutation_null_counts.csv'
    validate_saved_summaries(original.parent)
    original_hash = sha256(original)
    metadata = json.loads((ROOT/'results/diet_code/permutation_run.json').read_text())
    independent = seed is not None
    if independent and seed == metadata['seed']:
        raise ValueError('The independent run must use a different seed')
    selected_seed = seed if independent else metadata['seed']
    previous_root = permutation.ROOT
    try:
        with tempfile.TemporaryDirectory(prefix='vegandiet-audit-replay-') as temp:
            replay = Path(temp)
            (replay/'data').symlink_to(ROOT/'data',target_is_directory=True)
            out = replay/'results/diet_code'
            out.mkdir(parents=True)
            for name in ['before_limma.csv.gz','after_limma.csv.gz']:
                (out/name).symlink_to(ROOT/'results/diet_code'/name)
            permutation.ROOT = replay
            permutation.run('diet_code',metadata['repetitions'],selected_seed,
                            metadata['batch_size'],metadata['workers'])
            actual = pd.read_csv(out/'permutation_null_counts.csv')
            if not independent:
                np.testing.assert_array_equal(actual,pd.read_csv(original))
            observed = pd.read_csv(ROOT/'results/diet_code/limma_summary.csv').set_index('time')
            observed_change = int(observed.loc['after','highlighted']-observed.loc['before','highlighted'])
            extreme = int(((actual.after-actual.before)>=observed_change).sum())
            report = {
                'method':('independent Monte Carlo check; seed specified before fitting' if independent
                          else 'full replay of current production code in an isolated output directory'),
                'repetitions':len(actual),'seed':selected_seed,
                'original_count_sha256':original_hash,
                'observed_change':observed_change,'shuffles_ge_observed_change':extreme,
                'p_upper':(extreme+1)/(len(actual)+1),
                'elapsed_seconds':round(time.monotonic()-started,2),
                'numpy_version':np.__version__,'rng':'numpy.default_rng (PCG64)',
                'permutation_script_sha256':sha256(ROOT/'scripts/permutation.py'),
            }
            if not independent:
                report['all_paired_counts_equal'] = True
    finally:
        permutation.ROOT = previous_root
    assert sha256(original) == original_hash
    audit_out = ROOT/'results/audit'
    audit_out.mkdir(exist_ok=True)
    if independent:
        actual.to_csv(audit_out/'independent_seed_counts.csv',index=False)
    name = 'independent_seed_summary.json' if independent else 'full_replay_summary.json'
    (audit_out/name).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


def main():
    started = time.monotonic()
    out = ROOT/'results/audit'
    out.mkdir(exist_ok=True)
    raw = ROOT/'data/raw'
    prepared = ROOT/'data/processed'
    original = ROOT/'results/diet_code'
    manifest = json.loads((raw/'manifest.json').read_text())
    source_hashes = {}
    for entry in manifest:
        source_hashes[entry['file']] = sha256(raw/entry['file'])
        assert source_hashes[entry['file']] == entry['sha256'], entry['file']
    print('Source manifest hashes verified.', flush=True)

    # Independent annotation parsing; match each matrix column by array ID.
    with (raw/'E-MTAB-17025.sdrf.txt').open() as handle:
        source_rows = list(csv.DictReader(handle, delimiter='\t'))
    rows_by_array = {}
    for row in source_rows:
        array = 'ID_'+row['Array Data File'].rsplit('_',1)[0]
        rows_by_array.setdefault(array, []).append(row)
    samples = pd.read_csv(prepared/'samples.csv', dtype=str)
    beta = np.memmap(prepared/'beta.bin',dtype='<f8',mode='r',shape=(812934,96))
    cpg = (prepared/'cpg_ids.txt').read_text().splitlines()
    digest = hashlib.sha256()
    with (raw/'b_values4.txt').open('rb') as handle:
        header = handle.readline()
        digest.update(header)
        arrays = [x.strip('"') for x in header.decode().strip().split('/t')]
        assert arrays == samples.array_id.tolist()
        assert set(arrays) == set(rows_by_array)
        mapping = {'sample_id':'Source Name', 'participant':'Characteristics[individual]',
                   'time':'Factor Value[time]', 'diet_code':'Characteristics[diet]',
                   'diet_named':'Characteristics[disease]', 'sex':'Characteristics[sex]'}
        for i,array in enumerate(arrays):
            assert len(rows_by_array[array]) == 2
            for source in rows_by_array[array]:
                for field,key in mapping.items():
                    assert samples.iloc[i][field] == source[key], (array,field)
                assert float(samples.iloc[i].age) == float(source['Characteristics[age]'])
                assert source['Characteristics[disease]'] == source['Factor Value[disease]']
        for i,line in enumerate(handle):
            digest.update(line)
            parts = line.decode().strip().split('/t')
            assert parts[0].strip('"') == cpg[i]
            # Separate parsing route from prepare_data.py's np.fromstring.
            values = np.array([float(x) for x in parts[1:]])
            np.testing.assert_array_equal(values,beta[i])
            assert np.isfinite(values).all() and ((values>0)&(values<1)).all()
            if (i+1)%200000 == 0:
                print(f'Compared {i+1:,} source rows with the binary matrix.',flush=True)
    assert i+1 == len(cpg) == 812934 and len(set(cpg)) == len(cpg)
    assert digest.hexdigest() == source_hashes['b_values4.txt']
    assert (prepared/'beta.bin').stat().st_size == 812934*96*8
    assert samples.groupby('participant').time.apply(lambda s: sorted(s)==['0','1']).all()
    assert samples.groupby('participant')[['diet_code','diet_named','sex','age']].nunique().eq(1).all().all()
    print('All source values, sample mappings and paired labels verified.',flush=True)

    baseline = samples[samples.time=='0'].set_index('participant').sort_index()
    sexes = pd.crosstab(baseline.diet_code,baseline.sex)
    sexes.to_csv(out/'public_sex_counts.csv')
    consistency = validate_saved_summaries(original)
    counts = pd.read_csv(original/'permutation_null_counts.csv')
    changes = counts.after-counts.before
    saved = pd.read_csv(original/'permutation_null_changes.csv')
    np.testing.assert_array_equal(saved[['before','after']],counts[['before','after']])
    np.testing.assert_array_equal(saved.change,changes)
    summary = pd.read_csv(original/'change_test_summary.csv').iloc[0]
    observed_change = int(summary.observed_after-summary.observed_before)
    extreme = int((changes>=observed_change).sum())
    assert extreme == summary.shuffles_ge_observed
    assert (1+extreme)/(1+len(changes)) == summary.change_p_upper

    # Include observed, extremes, threshold neighbours and spread-out assignments.
    indices = {0, len(counts)-1, int(changes.idxmin()), int(changes.idxmax()), 99, 999, 4999}
    for value in [observed_change-1,observed_change,observed_change+1]:
        indices.add(int(np.flatnonzero(changes==value)[0]))
    metadata = json.loads((original/'permutation_run.json').read_text())
    rng = np.random.default_rng(metadata['seed'])
    labels = baseline.diet_code.eq('0').to_numpy()
    assignments = [{'assignment':'observed',**baseline.diet_code.to_dict()}]
    for i in range(len(counts)):
        shuffled = rng.permutation(labels)
        assert int(shuffled.sum()) == 24
        if i in indices:
            assignments.append({'assignment':str(i+1),
                **dict(zip(baseline.index,np.where(shuffled,'0','1')))})
    assignment_path = out/'reference_assignments.csv'
    pd.DataFrame(assignments).to_csv(assignment_path,index=False)
    subprocess.run(['bash',str(ROOT/'scripts/run_r.sh'),str(ROOT/'scripts/audit_limma.R'),
                    str(assignment_path),str(out/'independent_limma_counts.csv')],cwd=ROOT,check=True)
    ref = pd.read_csv(out/'independent_limma_counts.csv',dtype={'assignment':str})
    observed_ref = ref[ref.assignment=='observed'].iloc[0]
    observed_counts = pd.read_csv(original/'limma_summary.csv').set_index('time')
    for row in ref.itertuples():
        expected = (observed_counts.loc[['before','after'],'highlighted'].to_numpy()
                    if row.assignment=='observed' else counts.iloc[int(row.assignment)-1].to_numpy())
        np.testing.assert_array_equal([row.before,row.after],expected)
    errors = {k:float(v) for k,v in observed_ref.items() if k.startswith('max_abs_')}
    assert max(errors.values()) < 1e-8
    report = {
        'source_sha256':source_hashes,
        'binary_sha256':sha256(prepared/'beta.bin'),
        'source_values_checked':812934*96,
        'all_source_values_and_cpg_order_identical':True,
        'source_annotations_match_prepared_samples':True,
        'participants':len(baseline),
        'public_numeric_sex_counts':{str(k):v for k,v in sexes.to_dict(orient='index').items()},
        'source_group_crosstab':pd.crosstab(baseline.diet_code,baseline.diet_named).to_dict(orient='index'),
        'independent_r_model':'uncentered M-values, intercept + diet + visit + diet:visit',
        'independent_real_data_fits':len(ref),
        'all_independent_r_counts_match':True,
        'observed_max_absolute_errors':errors,
        'observed_change':observed_change,
        'shuffles_ge_observed':extreme,'permutations':len(counts),
        'p_upper':float(summary.change_p_upper),
        'summary_consistency':consistency,
        'permutation_count_sha256':sha256(original/'permutation_null_counts.csv'),
        'elapsed_seconds':round(time.monotonic()-started,2),
    }
    (out/'audit_summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--replay', action='store_true', help='also replay all saved shuffles')
    group.add_argument('--replay-only', action='store_true', help='only replay the saved shuffles')
    group.add_argument('--independent-seed', type=int,
                       help='run a separate Monte Carlo check without replacing the main results')
    args = parser.parse_args()
    if args.independent_seed is not None:
        replay_permutations(seed=args.independent_seed)
    elif not args.replay_only:
        main()
    if args.replay or args.replay_only:
        replay_permutations()

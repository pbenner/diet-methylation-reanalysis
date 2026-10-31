"""Explicit paths for archived annotations and the documented visit correction."""
SCENARIOS = ('corrected_visits', 'deposited')


def result_base(root, scenario):
    if scenario not in SCENARIOS:
        raise ValueError(f'Unknown annotation scenario: {scenario}')
    return root/'results/corrected_visits' if scenario == 'corrected_visits' else root/'results'


def paths(root, field='diet_code', scenario='corrected_visits'):
    folder = result_base(root, scenario)/field
    name = 'samples_corrected_visits.csv' if scenario == 'corrected_visits' else 'samples.csv'
    return root/'data/processed'/name, folder


def add_argument(parser):
    parser.add_argument('--scenario', choices=SCENARIOS, default='corrected_visits')

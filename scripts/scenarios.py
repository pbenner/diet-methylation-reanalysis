"""Explicit paths for the current blocked analysis and archived analysis scenarios."""
SCENARIOS = ('blocked_randomization', 'corrected_visits', 'deposited')
DEFAULT_SCENARIO = 'blocked_randomization'


def result_base(root, scenario):
    if scenario not in SCENARIOS:
        raise ValueError(f'Unknown annotation scenario: {scenario}')
    return root/'results'/scenario if scenario != 'deposited' else root/'results'


def paths(root, field='diet_code', scenario=DEFAULT_SCENARIO):
    if scenario == 'blocked_randomization' and field != 'diet_code':
        raise ValueError('The current analysis uses only author-verified numeric diet codes')
    folder = result_base(root, scenario)/field
    name = 'samples_corrected_visits.csv' if scenario != 'deposited' else 'samples.csv'
    return root/'data/processed'/name, folder


def add_argument(parser):
    parser.add_argument('--scenario', choices=SCENARIOS, default=DEFAULT_SCENARIO)

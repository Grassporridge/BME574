"""Run every antibiotic analysis script into separate output folders."""
import argparse
from pathlib import Path
import subprocess
import sys

SCRIPTS = [
    ('control', 'run_control.py'),
    ('dose_response', 'run_dose_response.py'),
    ('decay', 'run_decay_analysis.py'),
    ('resistance_selection', 'run_resistance_selection.py'),
    ('resistance_cost', 'run_resistance_cost_analysis.py'),
    ('growth_death', 'run_growth_death_analysis.py'),
    ('resource_interaction', 'run_resource_interaction.py'),
    ('sensitivity', 'run_sensitivity_analysis.py'),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--output', default='antibiotic_analysis_outputs')
    parser.add_argument('--sensitivity-samples', type=int, default=40)
    args = parser.parse_args()

    root = Path(args.output)
    root.mkdir(parents=True, exist_ok=True)
    here = Path(__file__).resolve().parent

    for name, script in SCRIPTS:
        cmd = [sys.executable, str(here / script), '--seed', str(args.seed), '--output', str(root / name)]
        if name == 'sensitivity':
            cmd += ['--samples', str(args.sensitivity_samples)]
        print('\n' + '=' * 72)
        print(f'Running {name}: {" ".join(cmd)}')
        subprocess.run(cmd, check=True, cwd=str(here))

    print('\nAll analyses completed.')
    print(f'Results: {root.resolve()}')


if __name__ == '__main__':
    main()

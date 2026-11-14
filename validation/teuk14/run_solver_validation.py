#!/usr/bin/env python3
"""Validate the type-14 CPU solver without a type-13 reference binary."""

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile

from run_validation import change

ROOT = Path(__file__).resolve().parents[2]


def number(log, label):
    match = re.search(re.escape(label) + r'\s+([-+0-9.eE]+)', log)
    if not match:
        raise AssertionError(f'Missing diagnostic: {label}')
    return float(match[1])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=Path, default=ROOT/'build/BSSN_GR/bssnSolver')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--ranks', type=int, default=2)
    parser.add_argument('--timeout', type=int, default=600)
    args = parser.parse_args()
    if args.ranks < 1 or args.timeout < 1:
        parser.error('--ranks and --timeout must be positive')
    output = (args.output or Path(tempfile.mkdtemp(prefix='teuk14-solver-'))).resolve()
    output.mkdir(parents=True, exist_ok=True)
    template = Path(__file__).with_name('uniform.toml').read_text()
    cases = {
        'uniform': ({}, args.ranks),
        'flat': ({'AMP': '0.0'}, args.ranks),
        'inactive': ({'BSSN_DENDRO_GRAIN_SZ': '100000'}, max(4, args.ranks)),
        'forced_failure': ({'TEUK_HAM_MAX_ITER': '1'}, args.ranks),
        'roundoff_guard': ({'TEUK_HAM_TOL': '1e-16', 'TEUK_HAM_MAX_ITER': '500'}, args.ranks),
    }
    summary = {}
    environment = dict(os.environ)
    environment.setdefault('OMP_NUM_THREADS', '1')
    for name, (overrides, ranks) in cases.items():
        run = output/name
        run.mkdir(exist_ok=True)
        for directory in ('cp', 'vtu', 'prof'):
            (run/directory).mkdir(exist_ok=True)
        text = template
        for key, value in overrides.items():
            text = change(text, key, value)
        (run/'input.toml').write_text(text)
        command = ['mpirun', '-np', str(ranks), str(args.solver.resolve()), 'input.toml', '1']
        with (run/'run.log').open('w') as log_file:
            result = subprocess.run(command, cwd=run, env=environment,
                                    stdout=log_file, stderr=subprocess.STDOUT,
                                    timeout=args.timeout)
        log = (run/'run.log').read_text(errors='replace')
        final = next((line for line in log.splitlines()
                      if line.startswith('TEUK_HAM final reason:')), '')
        accepted = '; accepted 1' in final
        true_relative = number(final, 'true relative residual')
        if name == 'forced_failure':
            assert result.returncode != 0 and not accepted and 'iteration limit' in final, log
        elif name == 'roundoff_guard':
            # Roundoff need not prevent convergence on every platform. It must
            # never produce acceptance above the requested true tolerance.
            assert (result.returncode == 0) == accepted, log
            assert not accepted or true_relative <= 1e-16, final
        else:
            assert result.returncode == 0 and accepted, log
            assert true_relative <= 1e-10, final
            threshold = 1e-11 if name == 'flat' else 3e-7
            assert number(log, 'TYPE14 solved C_HAM L2') < threshold, log
            assert number(log, 'TYPE14 min/max psi') > 0.0, log
            assert number(log, 'TYPE14 Gamma constraint L2/max') < 1e-12, log
            for axis in range(3):
                assert number(log, f'TYPE14 C_MOM{axis} L2') < 1e-12, log
            if name == 'inactive':
                assert 0 < number(log, 'TYPE14 active MPI ranks') < ranks, log
                reference = summary['uniform']['hamiltonian_rms']
                assert abs(number(log, 'TYPE14 solved C_HAM L2')/reference - 1) < 1e-5, log
        summary[name] = {'returncode': result.returncode, 'diagnostic': final}
        if name in ('uniform', 'flat', 'inactive'):
            summary[name]['hamiltonian_rms'] = number(log, 'TYPE14 solved C_HAM L2')
        (output/'results.json').write_text(json.dumps(summary, indent=2) + '\n')
        print(f'{name}: passed', flush=True)
    print(f'All solver checks passed; results: {output}')


if __name__ == '__main__':
    main()

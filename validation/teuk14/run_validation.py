#!/usr/bin/env python3
"""Fixed-mesh numerical regression; run after building bssnSolver.

Writes generated parameters and logs only under --output. The reference binary
is read without modification. Run from any working directory.
"""
import argparse
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[2]

def change(text, key, value):
    pattern = rf'^{re.escape(key)}\s*=.*$'
    if not re.search(pattern, text, re.M):
        raise ValueError(f'Missing template parameter {key}')
    return re.sub(pattern, lambda _: f'{key} = {value}', text, flags=re.M)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=Path, default=ROOT/'build/BSSN_GR/bssnSolver')
    parser.add_argument('--reference', type=Path, default=ROOT.parent/'teukolsky_solved_id.bin')
    parser.add_argument('--output', type=Path, default=Path('/tmp/dendro-teuk14-validation'))
    parser.add_argument('--ranks', type=int, default=2)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    args.output = args.output.resolve()
    template = (Path(__file__).parent/'uniform.toml').read_text()
    template = change(template, 'TEUK_SOLVED_ID_FILE', json.dumps(str(args.reference.resolve())))
    cases = {
        '49': ({}, args.ranks, True),
        '97': ({'BSSN_MAXDEPTH': '6'}, args.ranks, True),
        'flat': ({'AMP': '0.0'}, args.ranks, False),
        'inactive': ({'BSSN_DENDRO_GRAIN_SZ': '100000'}, max(4,args.ranks), False),
        'failure': ({'TEUK_HAM_MAX_ITER': '1'}, args.ranks, False),
    }
    logs = {}
    for name,(overrides,ranks,compare) in cases.items():
        text = template
        for key,value in overrides.items():
            text = change(text,key,value)
        par = args.output/f'{name}.toml'
        par.write_text(text)
        command = ['mpirun','-np',str(ranks),str(args.solver.resolve()),str(par),'1']
        if compare:
            command += ['--teuk-compare']
        with (args.output/f'{name}.log').open('w') as log:
            result = subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,timeout=600)
        logs[name] = (args.output/f'{name}.log').read_text()
        if name == 'failure':
            assert result.returncode != 0 and 'failed true-residual/positivity checks' in logs[name]
        else:
            assert result.returncode == 0, f'{name} failed; see its log'
        print(name, 'passed', flush=True)
    def number(case, label):
        match = re.search(re.escape(label)+r'\s+([-+0-9.eE]+)',logs[case])
        assert match, (case,label)
        return float(match[1])
    coarse = number('49','TYPE14 solved C_HAM L2')
    fine = number('97','TYPE14 solved C_HAM L2')
    assert coarse/fine>20, 'Hamiltonian does not improve with refinement'
    assert number('flat','TYPE14 solved C_HAM L2')<1e-11
    assert 0 < number('inactive','TYPE14 active MPI ranks') < max(4,args.ranks)
    assert abs(number('inactive','TYPE14 solved C_HAM L2')/coarse-1)<1e-5
    for case in ('49','97','flat','inactive'):
        assert number(case,'TYPE14 elliptic residual L2')<1e-11
        assert number(case,'TYPE14 Gamma constraint L2/max')<1e-12
        for i in range(3):
            assert number(case,f'TYPE14 C_MOM{i} L2')<1e-12
    assert number('97','COMPARE core type 13 C_HAM L2/max') > 1000*fine
    # Matching reference nodes must agree in the interior truncation scale.
    assert abs(number('49','COMPARE core type 13 C_HAM L2/max')/
               number('49','COMPARE core type 14 C_HAM L2/max')-1)<0.02
    summary = '\n'.join(f'[{case}]\n'+'\n'.join(line for line in log.splitlines()
        if line.startswith(('TYPE14','COMPARE','Type 14 failed'))) for case,log in logs.items())
    (args.output/'summary.txt').write_text(summary+'\n')
    print(f'Numerical checks passed; results: {args.output}')

if __name__ == '__main__':
    main()

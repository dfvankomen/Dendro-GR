# Porting testing scripts to HPC sites

Site-specific gotchas for the runners under `scripts/` (e.g. `hybrid_scaling/`).
When a script fails on a new machine it's almost always one of: modules, a
hardcoded path, or the MPI fabric. Patterns below; per-site recipes at the bottom.

## The three portability rules

1. **Modules** — don't cram a full module sequence into one `MODULES=...` var; a
   `module use <path>` step (custom module trees) can't fit `module load $VAR`, and
   names/versions differ per site (`intel-oneapi-mkl` vs `mkl`, `openmpi/5.0.3` vs
   `5.0.9`). Write the sequence explicitly per site and give a `MODULES=skip`
   escape hatch for users who load modules in an sbatch wrapper or `~/.bashrc`.

2. **Paths** — no personal defaults (`$HOME/research/...` only exists on one
   machine). Derive from the repo where possible, otherwise leave the default empty
   and let the preflight check fail loudly (`ERROR: set DENDROLIB_DIR ...`). A
   wrong-but-present default is worse than a clear failure.

3. **MPI fabric** — the wrong provider silently falls back to TCP and inflates
   inter-node comm times (the very thing a scaling benchmark measures). Set the
   fabric vars per site and make them overridable (`FI_PROVIDER="${FI_PROVIDER:-opx}"`).

## Site recipes

### Stampede3 (TACC) — Omni-Path fabric

```bash
module purge
module load autotools/1.4 cmake/4.1.1 python/3.12.11 gcc/15.1.0 mkl/25.1 gsl/2.8
module use /scratch/projects/compilers/modulefiles     # openmpi lives off the default tree
module load openmpi/5.0.9
export OMPI_MCA_mtl=ofi FI_PROVIDER=opx                 # REQUIRED: select Omni-Path, not TCP
# dendrolib is not under $HOME/research here, e.g.:
export DENDROLIB_DIR=$HOME/Projects/dendro/DVK_experimental/Dendro-5.01
```

Then run with `MODULES=skip` so the script doesn't re-`module load` its own stack:
```bash
MODULES=skip CPU_ARCH=sapphirerapids CORES_PER_NODE=112 \
  sbatch -A MYACCT -p spr --nodes=4 hybrid_scaling/run_hybrid_scaling.sh
```

### CHPC (Utah) — csl / rom nodes

The `hybrid_scaling/` defaults target this: `MODULES="gcc/15.1.0
intel-oneapi-mkl/... openmpi/..."`, `DENDROLIB_DIR=$HOME/research/dendrolib_dfvk_copy`,
`CPU_ARCH=cascadelake` (csl) or `znver2` (rom). See `hybrid_scaling/README.md`.

## Summary

| symptom | cause | fix |
|---------|-------|-----|
| `module load` fails | needs `module use <path>` first; wrong names/versions | inline the site sequence; `MODULES=skip` escape hatch |
| path / dir not found | personal default (`$HOME/research/...`) | repo-relative, or empty default + loud preflight failure |
| inter-node comm looks slow | MPI fell back to TCP | set fabric vars (Stampede3: `OMPI_MCA_mtl=ofi FI_PROVIDER=opx`), make overridable |

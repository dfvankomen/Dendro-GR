#!/usr/bin/env bash
# Pure-MPI vs hybrid MPI+OpenMP scaling for DendroGR BSSN. See README.md.
# Builds bssnScalingBench once, then sweeps OpenMP threads/rank at fixed total
# cores (T=1 is pure-MPI); parses the per-step JSONL into a CSV. Run multi-node.

#SBATCH --job-name=bssn_hybrid_scaling
#SBATCH --nodes=2
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --output=hybrid_scaling_%j.out
# account/partition/constraint differ per allocation -> pass on the sbatch line.

set -uo pipefail

# --- config (all overridable via env) ----------------------------------------
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
REPO="$HERE"; while [[ "$REPO" != / && ! -f "$REPO/CMakeLists.txt" ]]; do REPO="$(dirname "$REPO")"; done

# Site preset: modules + arch + cores/node + dendrolib path + MPI fabric. Default
# stampede3 (TACC); 'chpc' for Utah csl nodes; 'none' if you load modules yourself.
# Any individual value below is still overridable by its own env var. See PORTING.md.
SITE="${SITE:-stampede3}"
case "$SITE" in
  stampede3)                                       # Sapphire Rapids, Omni-Path fabric
    : "${CPU_ARCH:=sapphirerapids}"; SITE_CORES=112
    : "${DENDROLIB_DIR:=$HOME/Projects/dendro/DVK_experimental/Dendro-5.01}"
    export OMPI_MCA_mtl="${OMPI_MCA_mtl:-ofi}" FI_PROVIDER="${FI_PROVIDER:-opx}" ;;  # select Omni-Path, not TCP
  chpc)                                            # Utah CHPC, csl (Cascade Lake)
    : "${CPU_ARCH:=cascadelake}"; SITE_CORES=40
    : "${DENDROLIB_DIR:=$HOME/research/dendrolib_dfvk_copy}" ;;
  none|skip) SITE_CORES="" ;;                      # no presets; load modules yourself
  *) echo "ERROR: unknown SITE=$SITE (stampede3|chpc|none)"; exit 1 ;;
esac
CPU_ARCH="${CPU_ARCH:-native}"
DENDROLIB_DIR="${DENDROLIB_DIR:-}"                 # from site preset or env; must have the hybrid path (checked below)
CORES_PER_NODE="${CORES_PER_NODE:-${SLURM_CPUS_ON_NODE:-${SITE_CORES:-$(nproc)}}}"

BUILD_ROOT="${BUILD_ROOT:-$HERE}"                  # put on shared scratch for multi-node
BUILD_DIR="${BUILD_DIR:-$BUILD_ROOT/build_hybrid}"
CASCADE_FLAG="${CASCADE_FLAG:-}"                   # optional, e.g. -DBSSN_USE_CASCADE_AVX512_FUSED=ON
JOBS="${JOBS:-$(nproc)}"
NNODES="${NNODES:-${SLURM_NNODES:-1}}"
TOTAL_CORES=$(( CORES_PER_NODE * NNODES ))
MPI_LAUNCH="${MPI_LAUNCH:-mpirun}"       # or srun (set SRUN_MPI plugin)
SRUN_MPI="${SRUN_MPI:-pmix}"

THREADS_LIST="${THREADS_LIST:-1 2 4}"    # threads/rank; 1 = pure-MPI. Keep each a per-socket divisor.

# Default = the self-contained BBH grid here. NOT BSSN_GR/pars/q1.par.toml: it needs
# an external TwoPunctures (tp_*) file and segfaults without one.
PARFILE="${PARFILE:-$HERE/q1.scaling.par.toml}"
if [[ "$PARFILE" != /* ]]; then
  [[ -f "$HERE/$PARFILE" ]] && PARFILE="$HERE/$PARFILE" || PARFILE="$REPO/BSSN_GR/pars/$PARFILE"
fi
GRID="${GRID:-bbh}"; LEV="${LEV:-7}"     # uniform grid blocks ~= 8^(LEV-3); use it to saturate many nodes
STEPS="${STEPS:-10}"; WARMUP="${WARMUP:-2}"

STAMP="${SLURM_JOB_ID:-local}"
OUTDIR="${OUTDIR:-$HERE/results_${STAMP}}"
OUT_CSV="${OUT_CSV:-$OUTDIR/hybrid_scaling_${STAMP}.csv}"

# --- preflight ---------------------------------------------------------------
mkdir -p "$OUTDIR"
if [[ "$SITE" != none && "$SITE" != skip ]] && command -v module >/dev/null 2>&1; then
  module purge 2>/dev/null || true
  case "$SITE" in
    stampede3)   # openmpi is off the default module tree -> needs 'module use' first
      module load autotools/1.4 cmake/4.1.1 python/3.12.11 gcc/15.1.0 mkl/25.1 gsl/2.8 \
        && module use /scratch/projects/compilers/modulefiles \
        && module load openmpi/5.0.9 ;;
    chpc)
      module load gcc/15.1.0 intel-oneapi-mkl/2025.3.1 openmpi/5.0.3 ;;
  esac || { echo "ERROR: module load failed for SITE=$SITE (try SITE=none if you load modules yourself)"; exit 1; }
fi
command -v cmake   >/dev/null || { echo "ERROR: cmake not found";   exit 1; }
command -v python3 >/dev/null || { echo "ERROR: python3 not found"; exit 1; }
[[ -f "$PARFILE" ]]     || { echo "ERROR: parfile not found: $PARFILE"; exit 1; }
[[ -d "$DENDROLIB_DIR" ]] || { echo "ERROR: DENDROLIB_DIR not set/found: '$DENDROLIB_DIR' (set DENDROLIB_DIR, or SITE with a preset)"; exit 1; }

echo "site=$SITE  nodes=$NNODES cores/node=$CORES_PER_NODE total=$TOTAL_CORES  launch=$MPI_LAUNCH arch=$CPU_ARCH"
echo "sweep T=[$THREADS_LIST]  parfile=$(basename "$PARFILE") grid=$GRID steps=$STEPS  -> $OUTDIR"

# --- build once (hybrid path + profile counters for the JSONL) ---------------
BIN="$BUILD_DIR/BSSN_GR/bssnScalingBench"
if [[ ! -x "$BIN" ]]; then
  echo "## building (once) ..."
  cmake -S "$REPO" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release -DCPU_ARCH="$CPU_ARCH" \
        -DDENDRO_dendrolib_DIR="$DENDROLIB_DIR" -DDENDRO_HYBRID_OMP=ON \
        -DENABLE_DENDRO_PROFILE_COUNTERS=ON $CASCADE_FLAG >"$OUTDIR/configure.log" 2>&1 \
    || { echo "## CONFIGURE FAILED -- $OUTDIR/configure.log"; tail -20 "$OUTDIR/configure.log"; exit 1; }
  cmake --build "$BUILD_DIR" --target bssnScalingBench -j"$JOBS" >"$OUTDIR/build.log" 2>&1 \
    || { echo "## BUILD FAILED -- $OUTDIR/build.log"; tail -20 "$OUTDIR/build.log"; exit 1; }
else
  echo "## reusing build $BIN (rm -rf $BUILD_DIR to rebuild)"
fi

# --- launch one split: R ranks x T threads. Pinning is what makes hybrid work. -
run_split() {
  local R="$1" T="$2" RPN="$3" PFX="$4"
  export OMP_NUM_THREADS="$T" OMP_PROC_BIND=close OMP_PLACES=cores
  local args=("$PARFILE" --grid "$GRID" --lev "$LEV" --steps "$STEPS" --warmup "$WARMUP" --prefix "$PFX")
  if [[ "$MPI_LAUNCH" == srun ]]; then
    srun --mpi="$SRUN_MPI" --nodes="$NNODES" --ntasks="$R" --ntasks-per-node="$RPN" \
         --cpus-per-task="$T" --cpu-bind=cores --distribution=block:block "$BIN" "${args[@]}"
  else
    mpirun --np "$R" --map-by "ppr:${RPN}:node:pe=${T}" --bind-to core "$BIN" "${args[@]}"
  fi
}

# --- sweep -------------------------------------------------------------------
for T in $THREADS_LIST; do
  (( CORES_PER_NODE % T == 0 )) || { echo "## skip T=$T (does not divide $CORES_PER_NODE cores/node)"; continue; }
  RPN=$(( CORES_PER_NODE / T )); R=$(( RPN * NNODES ))
  PFX="$OUTDIR/run_N${NNODES}_r${R}_t${T}"; rm -f "${PFX}_steps.jsonl"
  echo "===== T=$T ($([[ $T == 1 ]] && echo pure-MPI || echo hybrid)) ranks=$R ranks/node=$RPN ====="
  run_split "$R" "$T" "$RPN" "$PFX" >"${PFX}.log" 2>&1
  [[ -s "${PFX}_steps.jsonl" ]] && echo "## ok" || { echo "## FAILED (no JSONL):"; tail -15 "${PFX}.log"; }
done

# --- parse -> CSV ------------------------------------------------------------
echo ""; echo "## parsing -> $OUT_CSV"
python3 - "$OUTDIR" "$OUT_CSV" <<'PYEOF'
import glob, json, os, re, sys
outdir, out_csv = sys.argv[1], sys.argv[2]
PHASES = ["rk_step", "rhs", "deriv", "unzip", "unzip_wcomm", "zip"]
mean = lambda xs: sum(xs) / len(xs) if xs else float("nan")

rows = []
for path in sorted(glob.glob(os.path.join(outdir, "run_N*_steps.jsonl"))):
    m = re.search(r"run_N(\d+)_r(\d+)_t(\d+)_steps\.jsonl$", os.path.basename(path))
    if not m: continue
    nodes, ranks, threads = map(int, m.groups())
    # per phase, mean over timed steps of the per-step max-across-ranks (the critical path)
    acc = {p: [] for p in PHASES}; nblocks = None
    for line in open(path):
        line = line.strip()
        if not line: continue
        try: rec = json.loads(line)
        except json.JSONDecodeError: continue
        ph = rec.get("phase", {})
        for p in PHASES:
            if p in ph and "max" in ph[p]: acc[p].append(float(ph[p]["max"]))
        if nblocks is None: nblocks = rec.get("num_blocks")
    v = {p: mean(acc[p]) for p in PHASES}
    rows.append(dict(mode="pure-MPI" if threads == 1 else "hybrid", nodes=nodes, ranks=ranks,
                     threads=threads, ranks_per_node=ranks // nodes if nodes else 0,
                     ghost_comm=v["unzip_wcomm"] - v["unzip"], num_blocks=nblocks,
                     blocks_per_rank=(nblocks / ranks) if (nblocks and ranks) else float("nan"), **v))

base = {r["nodes"]: r["rk_step"] for r in rows if r["threads"] == 1}
for r in rows:
    b = base.get(r["nodes"]); r["speedup"] = (b / r["rk_step"]) if (b and r["rk_step"] > 0) else float("nan")
rows.sort(key=lambda r: (r["nodes"], r["threads"]))

cols = ["mode", "nodes", "ranks", "threads", "ranks_per_node", "rk_step", "rhs", "deriv",
        "unzip", "unzip_wcomm", "ghost_comm", "zip", "num_blocks", "blocks_per_rank", "speedup"]
hdr = {"rk_step": "rk_step_s", "rhs": "rhs_s", "deriv": "deriv_s", "unzip": "unzip_s",
       "unzip_wcomm": "unzip_wcomm_s", "ghost_comm": "ghost_comm_s", "zip": "zip_s",
       "speedup": "speedup_vs_purempi"}
fmt = lambda x: (f"{x:.4f}" if x == x else "") if isinstance(x, float) else ("" if x is None else str(x))
with open(out_csv, "w") as fh:
    fh.write(",".join(hdr.get(c, c) for c in cols) + "\n")
    for r in rows: fh.write(",".join(fmt(r[c]) for c in cols) + "\n")
print(f"wrote {len(rows)} rows -> {out_csv}")

# saturation check: too few blocks/rank => imbalance dominates, results are noise
bprs = [r["blocks_per_rank"] for r in rows if r["blocks_per_rank"] == r["blocks_per_rank"]]
if bprs and min(bprs) < 4:
    print(f"\n  !! UNDER-SATURATED: as few as {min(bprs):.1f} blocks/rank -- results unreliable."
          f"\n     Use a bigger mesh: GRID=uniform LEV=7 (~4096 blocks) or LEV=8, or fewer nodes.")
elif bprs and min(bprs) < 16:
    print(f"\n  ! marginal: as few as {min(bprs):.1f} blocks/rank (prefer >=16). Try a larger GRID=uniform LEV.")
PYEOF

echo ""; echo "== $OUT_CSV =="
column -t -s, "$OUT_CSV" 2>/dev/null || cat "$OUT_CSV"
echo "(headline rk_step_s; ghost_comm_s = the inter-node comm hybrid shrinks; speedup>1 = hybrid wins."
echo " phase times carry profiler-barrier overhead -> read them relative, not absolute.)"

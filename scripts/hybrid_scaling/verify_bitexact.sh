#!/usr/bin/env bash
# Bit-exactness gate for mesh-construction changes (Phase B: threading
# buildE2EMap / buildE2NMap / performBlocksSetup / OCT2BLK).
#
# The invariant: the mesh is a deterministic function of the octree, so at a
# FIXED rank count its structure must not depend on the thread count. Any
# divergence is a race or a reordering -- and a reordered E2N map silently
# changes every numeric result forever, which is the failure this exists to stop.
#
# WHY RANKS ARE HELD FIXED: the mesh is legitimately rank-count dependent. The
# SFC partition gives each rank different elements, node numbering is rank-local,
# and getAllElements() includes ghosts. So comparing across a fixed-CORE sweep
# (R = cores/T) would show differences that are entirely correct. We vary ONLY
# OMP_NUM_THREADS, oversubscribing if need be -- this is a correctness gate, so
# timing does not matter.
#
# Digests are per construction stage, so a failure localizes: e2n differing while
# e2e matches points straight at buildE2NMap.
#
#   ./verify_bitexact.sh                  # default: R=4, T=1,2,4, depth 9
#   RANKS=8 THREADS_LIST="1 4" ./verify_bitexact.sh
#   BIN=<other build>/bssnScalingBench ./verify_bitexact.sh   # e.g. flag OFF vs ON
set -uo pipefail

REPO="${REPO:-/home/denv/research/dendrogr_dfvk}"
BIN="${BIN:-$REPO/build_hybrid/BSSN_GR/bssnScalingBench}"
BASE_PAR="${BASE_PAR:-$REPO/scripts/hybrid_scaling/q1.scaling.par.toml}"
OUTDIR="${OUTDIR:-$PWD/bitexact}"
RANKS="${RANKS:-4}"                       # FIXED across the sweep -- see above
THREADS_LIST="${THREADS_LIST:-1 2 4}"
DEPTH="${DEPTH:-9}"                       # small: correctness, not performance
STEPS="${STEPS:-2}"
WARMUP="${WARMUP:-1}"

mkdir -p "$OUTDIR"
[[ -x "$BIN" ]] || { echo "no binary at $BIN (build target bssnScalingBench)"; exit 2; }

PAR="$OUTDIR/verify.d${DEPTH}.par.toml"
sed "s/^BSSN_MAXDEPTH *=.*/BSSN_MAXDEPTH = $DEPTH/" "$BASE_PAR" > "$PAR"

echo "bit-exactness gate: R=$RANKS (fixed)  T=[$THREADS_LIST]  depth=$DEPTH  bin=$BIN"

ref=""; rc=0
for T in $THREADS_LIST; do
  f="$OUTDIR/fp_t${T}.txt"
  OMP_NUM_THREADS="$T" OMP_PROC_BIND=close OMP_PLACES=cores \
    mpirun --np "$RANKS" --oversubscribe "$BIN" "$PAR" \
      --grid bbh --steps "$STEPS" --warmup "$WARMUP" --fingerprint \
      --prefix "$OUTDIR/be_t${T}" > "$OUTDIR/run_t${T}.log" 2>&1
  grep "^\[fingerprint\]" "$OUTDIR/run_t${T}.log" | awk '{print $2,$3,$4}' > "$f"
  if [[ ! -s "$f" ]]; then
    echo "  T=$T: FAILED to produce digests (see $OUTDIR/run_t${T}.log)"; rc=1; continue
  fi
  # Guard against a digest that means "nothing was hashed" -- otherwise an
  # empty hash would sail through as a pass.
  if grep -q "cbf29ce484222325" "$f"; then
    echo "  T=$T: SUSPECT -- a digest is the un-hashed FNV offset (hashed nothing)"; rc=1
  fi
  if [[ -z "$ref" ]]; then
    ref="$f"; echo "  T=$T: reference"; sed 's/^/      /' "$f"
  elif diff -q "$ref" "$f" >/dev/null; then
    echo "  T=$T: MATCH"
  else
    echo "  T=$T: *** MISMATCH ***  (the differing stage names the culprit)"
    diff "$ref" "$f" | sed 's/^/      /'
    rc=1
  fi
done

echo
if (( rc == 0 )); then
  echo "PASS - mesh + evolved state are bit-identical across threads at R=$RANKS"
else
  echo "FAIL - see above. A stage listed here is not thread-invariant."
fi
exit $rc

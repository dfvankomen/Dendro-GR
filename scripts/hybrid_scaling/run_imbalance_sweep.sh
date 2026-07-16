#!/usr/bin/env bash
# Mechanism test: is the hybrid win really dynamic-schedule load balancing?
#
# Sweeps BSSN_MAXDEPTH (=> block count) x T at FIXED total cores, and records
# rhs load imbalance alongside speedup. The causal chain we are testing:
#
#   raising T pools T ranks' blocks under one schedule(dynamic,1) loop
#     -> rhs imbalance (max/mean across ranks) falls
#     -> rk_step falls
#   and that effect shrinks as blocks_per_rank rises, because a finer partition
#   has less quantization error left to recover.
#
# PRIMARY evidence is speedup_vs_T1 as a function of blocks_per_rank: it should
# be large at low bpr and collapse toward (or below) 1.0 as bpr grows. That is
# pure timing, with no sampling artifact.
#
# Imbalance is DIAGNOSTIC ONLY, because R and T are coupled here (R = CORES/T):
# as T rises the max is taken over FEWER ranks, so max/mean shrinks mechanically
# even with an identical work distribution -- at R=1 it is 1.000 by definition.
# So a fall in imbalance with T is NOT clean evidence for the mechanism. Read it
# across depths at FIXED T (where R is constant) instead; that comparison is
# sound and is what shows quantization imbalance decaying with bpr.
#
# Imbalance is read from rhs, NOT rk_step: the ghost exchange synchronizes
# ranks, so a rank that finishes early just blocks longer in comm and rk_step
# max/mean reads ~1.000 regardless of skew.
#
# Single-socket laptop caveat: absolute speedups are compressed by the 2-channel
# memory wall (see Single Node Split Sweep, ~4.5% at best). We are testing the
# SHAPE of imbalance-vs-T, not the magnitude of the speedup.
set -euo pipefail

REPO="${REPO:-/home/denv/research/dendrogr_dfvk}"
DENDROLIB="${DENDROLIB:-/home/denv/research/dendrolib_dfvk_copy}"
BUILD="${BUILD:-$REPO/build_hybrid}"
BIN="$BUILD/BSSN_GR/bssnScalingBench"
BASE_PAR="${BASE_PAR:-$REPO/scripts/hybrid_scaling/q1.scaling.par.toml}"
OUTDIR="${OUTDIR:-$PWD/imbalance_sweep}"

CORES="${CORES:-8}"                        # physical cores
DEPTHS="${DEPTHS:-8 9 10 11 12 13 14}"     # -> block count -> blocks_per_rank
THREADS_LIST="${THREADS_LIST:-1 2 4}"
STEPS="${STEPS:-6}"
WARMUP="${WARMUP:-2}"

mkdir -p "$OUTDIR"
echo "sweep: cores=$CORES depths=[$DEPTHS] threads=[$THREADS_LIST] steps=$STEPS -> $OUTDIR"

[[ -x "$BIN" ]] || {
  echo "building $BIN ..."
  cmake -S "$REPO" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release -DCPU_ARCH=native \
        -DDENDRO_dendrolib_DIR="$DENDROLIB" -DDENDRO_HYBRID_OMP=ON \
        -DENABLE_DENDRO_PROFILE_COUNTERS=ON -DOCT2BLK_COARSEST_LEV=0 >/dev/null
  cmake --build "$BUILD" --target bssnScalingBench -j"$(nproc)" >/dev/null
}

for D in $DEPTHS; do
  PAR="$OUTDIR/q1.d${D}.par.toml"
  sed "s/^BSSN_MAXDEPTH *=.*/BSSN_MAXDEPTH = $D/" "$BASE_PAR" > "$PAR"
  for T in $THREADS_LIST; do
    R=$(( CORES / T ))
    # R=1 makes max/mean across ranks trivially 1.000 -- no information.
    (( R >= 2 )) || { echo "  depth=$D  T=$T  R=$R -- skipped (R<2, imbalance degenerate)"; continue; }
    PFX="$OUTDIR/d${D}_t${T}"
    # blocks_per_rank is what we are sweeping, so record the bench's own count.
    LOG="${PFX}.log"
    echo "  depth=$D  T=$T  R=$R"
    OMP_NUM_THREADS="$T" OMP_PROC_BIND=close OMP_PLACES=cores \
      mpirun --np "$R" --map-by "ppr:${R}:node:pe=${T}" --bind-to core \
        "$BIN" "$PAR" --grid bbh --steps "$STEPS" --warmup "$WARMUP" \
        --prefix "$PFX" > "$LOG" 2>&1 || { echo "    FAILED (see $LOG)"; continue; }
    grep -m1 "\[bench\] elements=" "$LOG" | sed 's/^/    /' || true
  done
done

python3 - "$OUTDIR" <<'PY'
import glob, json, os, re, sys
outdir = sys.argv[1]
mean = lambda xs: sum(xs)/len(xs) if xs else float("nan")
rows = []
for path in sorted(glob.glob(os.path.join(outdir, "d*_t*_steps.jsonl"))):
    m = re.search(r"d(\d+)_t(\d+)_steps\.jsonl$", os.path.basename(path))
    if not m: continue
    depth, threads = int(m.group(1)), int(m.group(2))
    rk, rhs, imb, nb, ranks, omp = [], [], [], None, None, None
    for line in open(path):
        line = line.strip()
        if not line: continue
        try: rec = json.loads(line)
        except json.JSONDecodeError: continue
        ph = rec.get("phase", {})
        if "rk_step" in ph: rk.append(float(ph["rk_step"]["max"]))
        r = ph.get("rhs")
        if r:
            rhs.append(float(r["max"]))
            if r.get("mean", 0) > 0: imb.append(float(r["max"])/float(r["mean"]))
        if nb is None: nb = rec.get("num_blocks")
        # active_npes, not CORES/T: the mesh may not span every rank
        if ranks is None: ranks = rec.get("active_npes")
        if omp is None: omp = rec.get("omp_threads")
    if not rk:
        print(f"  !! no usable records in {os.path.basename(path)} -- skipped")
        continue
    # CANARY: omp_threads IS BSSN_HYBRID_NTHREADS. If it doesn't track
    # OMP_NUM_THREADS the BSSN omp regions ran num_threads(1) and this row is
    # meaningless -- it measures "more blocks per rank, same one thread". That
    # failure is bit-exact and silent (see 9bb8f45), so nothing else catches it.
    if omp is not None and omp != threads:
        print(f"  !! BOGUS: depth={depth} T={threads} but omp_threads={omp} "
              f"-- hybrid path DEAD, row discarded. Check BSSN_HYBRID_NTHREADS.")
        continue
    rows.append(dict(depth=depth, threads=threads, ranks=ranks, blocks=nb,
                     bpr=(nb/ranks if (nb and ranks) else float("nan")),
                     rk=mean(rk), rhs=mean(rhs), imb=mean(imb)))

# baseline per depth = T=1; bpr quoted at T=1 (the partition being fixed)
base   = {r["depth"]: r["rk"]  for r in rows if r["threads"] == 1}
bpr_t1 = {r["depth"]: r["bpr"] for r in rows if r["threads"] == 1}
for r in rows:
    b = base.get(r["depth"])
    r["speedup"]  = (b / r["rk"]) if (b and r["rk"] > 0) else float("nan")
    r["bpr_t1"]   = bpr_t1.get(r["depth"], float("nan"))
rows.sort(key=lambda r: (r["depth"], r["threads"]))

cols = ["depth","blocks","bpr_t1","ranks","threads","bpr","rk","rhs","imb","speedup"]
hdr  = {"bpr_t1":"blocks_per_rank_at_T1","bpr":"blocks_per_rank","rk":"rk_step_s",
        "rhs":"rhs_s","imb":"rhs_imbalance_max_over_mean","speedup":"speedup_vs_T1"}
fmt  = lambda x: (f"{x:.4f}" if x == x else "") if isinstance(x, float) else ("" if x is None else str(x))
csv = os.path.join(outdir, "imbalance_sweep.csv")
with open(csv, "w") as fh:
    fh.write(",".join(hdr.get(c, c) for c in cols) + "\n")
    for r in rows: fh.write(",".join(fmt(r[c]) for c in cols) + "\n")
print(f"\nwrote {len(rows)} rows -> {csv}\n")

# PRIMARY: does the hybrid speedup decay as blocks_per_rank rises?
print("PRIMARY -- speedup vs blocks_per_rank (should decay toward/below 1.0 as B rises):")
print(f"  {'B@T1':>8} {'T':>3} {'ranks':>6} {'speedup':>9} {'rk_step_s':>10}")
for d in sorted(set(r['depth'] for r in rows)):
    for r in [x for x in rows if x["depth"] == d]:
        print(f"  {r['bpr_t1']:8.1f} {r['threads']:3d} {r['ranks']:6d} "
              f"{r['speedup']:9.3f} {r['rk']:10.4f}")
print()
tops = sorted({r["threads"] for r in rows} - {1})
for T in tops:
    pts = sorted([r for r in rows if r["threads"] == T], key=lambda r: r["bpr_t1"])
    if len(pts) >= 2:
        print(f"  T={T}: speedup {pts[0]['speedup']:.3f} at B={pts[0]['bpr_t1']:.0f}"
              f"  ->  {pts[-1]['speedup']:.3f} at B={pts[-1]['bpr_t1']:.0f}"
              f"   ({'DECAYS as predicted' if pts[0]['speedup'] > pts[-1]['speedup'] else 'DOES NOT decay -- theory in trouble'})")
print()
# DIAGNOSTIC: compare across depths at FIXED T only (R constant => no sampling artifact).
print("DIAGNOSTIC -- rhs imbalance vs B at FIXED T (compare DOWN each block, never across T:")
print("             R=CORES/T, so max/mean shrinks with T mechanically):")
for T in sorted({r["threads"] for r in rows}):
    pts = sorted([r for r in rows if r["threads"] == T], key=lambda r: r["bpr_t1"])
    if not pts: continue
    print(f"  T={T} (R={pts[0]['ranks']}):  " +
          "  ".join(f"B={r['bpr_t1']:.0f}->{r['imb']:.3f}" for r in pts))
PY

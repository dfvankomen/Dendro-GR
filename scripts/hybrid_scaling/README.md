# Hybrid MPI + OpenMP scaling test (DendroGR / BSSN)

Does the hybrid **MPI + OpenMP** path beat **pure-MPI** when there's real inter-node
communication? `run_hybrid_scaling.sh` builds once, sweeps MPI×OpenMP splits at a
fixed core count, and writes a CSV. **Use as many nodes as you can** — single-node
doesn't exercise the interconnect, which is the point.

## Run it

```bash
cd scripts/hybrid_scaling
# SITE defaults to stampede3 (modules, arch, cores, dendrolib, MPI fabric all preset).
sbatch -A MYACCT -p spr --nodes=4 run_hybrid_scaling.sh            # Stampede3 (default)
SITE=chpc sbatch -A MYACCT -p soc-np -C csl --nodes=4 run_hybrid_scaling.sh   # CHPC csl
```

Output → `results_<jobid>/hybrid_scaling_<jobid>.csv` (plus per-run logs & build logs).

`SITE` sets the module stack, MPI-fabric vars, `CPU_ARCH`, `CORES_PER_NODE`, and
`DENDROLIB_DIR`. Use `SITE=none` if you load modules yourself. New site? Add a case
to the script and see [`../PORTING.md`](../PORTING.md).

## How it works

Total cores are held fixed (`cores/node × nodes`); the sweep varies OpenMP
threads/rank `T` (so ranks `R = totalcores/T`). `T=1` is pure-MPI; `T>1` is hybrid.
Same cores, same nodes → compute is constant and only the **number of MPI ranks
exchanging ghost cells** changes: bigger `T` = fewer ranks = fewer/larger inter-node
messages. Pure-MPI and hybrid are the *same binary* — pure-MPI is just
`OMP_NUM_THREADS=1`. The benchmark `bssnScalingBench` runs the full per-step compute
(no remesh/IO) and writes one JSON line per step.

## Reading the CSV

| column | meaning |
|--------|---------|
| `rk_step_s` | **headline** — mean wall time per RK step (lower = better) |
| `rhs_s` `deriv_s` `zip_s` | compute phases (should stay ~flat across the sweep) |
| `ghost_comm_s` | `unzip_wcomm − unzip` = the inter-node **communication** hybrid shrinks |
| `blocks_per_rank` | **saturation** — keep ≥ 16 or results are noise (script warns if low) |
| `speedup_vs_purempi` | pure-MPI rk_step / this rk_step (>1 → hybrid wins) |

**Hybrid winning** = flat compute, falling `ghost_comm_s`, rising speedup. If
`ghost_comm_s` is already tiny you're not comm-bound — use more nodes.

> Profiler barriers (needed to emit the JSONL) inflate the unzip phases, so read
> phase times as **relative** (pure-MPI vs hybrid), not absolute.

## Sizing the mesh

The mesh is fixed by the parfile and does **not** grow with node count, so the
busiest config (pure-MPI) needs `blocks_per_rank ≥ 16`. The default BBH grid is
~1044 blocks — good for ~1 node (the comm-bound regime where hybrid helps most), but
it **plateaus**: raising `BSSN_MAXDEPTH` past 14 adds no blocks. To saturate many
nodes, switch to a uniform grid (`blocks ≈ 8^(LEV−3)`):

`GRID=uniform` auto-builds with block fusion **off** (`OCT2BLK_LEV=31`,
one block per element) — otherwise a uniform grid collapses into a few giant fused
blocks and can't fill ranks; this is what makes `blocks ≈ 8^(LEV−3)` hold. `bbh`
keeps fusion **on** (`=0`), the production RHS layout.

No-fusion is memory-hungry: every one of the `8^LEV` blocks carries its own ghost
zone, so the unzipped arrays grow fast and OOM (as a `bad_alloc`, or a corrupt-pointer
segfault when it hits inside an OpenMP region). `MEM_CHECK` estimates per-node RAM and
aborts with a node count before that happens. Rough fit: `LEV=5` (~32K blocks) runs
2→8 nodes on one SKX/SPR node's worth of memory; `LEV=7` (~2M blocks) needs ~15+ nodes.

```bash
GRID=uniform LEV=7 sbatch --nodes=8 ... run_hybrid_scaling.sh   # ~4096 blocks (LEV=8 ~= 32768)
```

## Roofline mode

The compute phases (unzip, RHS) are memory-bandwidth-bound — they speed up as ranks/socket
drop because there's less bandwidth contention, not more parallelism. `PROFILE=roofline` tells
you how close you are to the memory wall, which decides whether *any* compute tuning is worth it:

```bash
module load likwid
PROFILE=roofline sbatch -A ACCT -p spr --nodes=1 run_hybrid_scaling.sh
```

It runs one config on a single node under `likwid-perfctr -g MEM_DP` (IMC counters see all node
memory traffic) and a `likwid-bench` STREAM run for the ceiling, then reports achieved GB/s, the
STREAM ceiling, **% of the bandwidth wall**, DP GFLOP/s, and operational intensity →
`roofline_<jobid>.csv`. Near the wall (~80%+) means only *fewer bytes* helps (ghost compression);
well below it means there's locality/prefetch headroom. Sweep `ROOFLINE_T` to see utilization vs split.

## Knobs (env vars)

| var | default | notes |
|-----|---------|-------|
| `SITE` | `stampede3` | preset bundle: modules + `CPU_ARCH` + `CORES_PER_NODE` + `DENDROLIB_DIR` + MPI fabric. Also `chpc`, or `none` (load modules yourself). |
| `THREADS_LIST` | `1 2 4` | threads/rank; `1` = pure-MPI. Keep each a divisor of the per-**socket** core count. |
| `GRID`/`LEV` | `bbh`/`7` | `uniform` + `LEV` to saturate many nodes (see above) |
| `OCT2BLK_LEV` | `31` uniform / `0` bbh | block fusion: `0` = fused big blocks (production RHS), `31` = no fusion, many small blocks. Baked into the build; changing it triggers a separate `build_hybrid_ob<N>/`. |
| `MEM_CHECK` | `on` | uniform-grid RAM preflight; aborts with a node count before an OOM. Set `off` to bypass the estimate. |
| `PROFILE` | *(none)* | `roofline` runs one single-node likwid measurement instead of the sweep (see below). |
| `ROOFLINE_T` | `1` | threads/rank for the roofline run; sweep it to see how bandwidth utilization changes with the split. |
| `STEPS`/`WARMUP` | `10`/`2` | timed / discarded RK steps |
| `MPI_LAUNCH` | `mpirun` | or `srun` (set `SRUN_MPI` plugin) |
| `CPU_ARCH` / `CORES_PER_NODE` / `DENDROLIB_DIR` | from `SITE` | override any of these to deviate from the site preset |
| `BUILD_ROOT` | this dir | put on **shared** scratch for multi-node so every node sees the binary |
| `CASCADE_FLAG` | *(none)* | optional RHS kernel, e.g. `-DBSSN_USE_CASCADE_AVX512_FUSED=ON` |

Don't point `PARFILE` at `BSSN_GR/pars/q1.par.toml` — it needs an external TwoPunctures
(`tp_*`) file and segfaults without one; the shipped `q1.scaling.par.toml` is
self-contained. Pinning is what makes hybrid pay off; if a split looks oddly slow,
check `mpirun --report-bindings`.

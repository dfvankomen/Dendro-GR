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
threads/rank `T` (so ranks `R = totalcores/T`). Same cores, same nodes → compute is
constant and only the **number of MPI ranks exchanging ghost cells** changes: bigger
`T` = fewer ranks = fewer/larger inter-node messages. The benchmark
`bssnScalingBench` runs the full per-step compute (no remesh/IO) and writes one JSON
line per step.

> **`T=1` is NOT pure MPI.** It is *hybrid-with-one-thread*, and the CSV labels it
> `hybrid-T1`. The binary is built `-DDENDRO_HYBRID_OMP=ON`, which forces
> `DENDRO_UNZIP_OMP_EFFECTIVE=ON` (`dendrolib/CMakeLists.txt:295`), so even at `T=1`
> it takes the OpenMP unzip path and pays that path's `all_dg` cost
> (`mesh.tcc:11295`) with no threading to offset it — measured at **~13% on unzip,
> ~6% on `rk_step`** (8 cores, depth 9, 2026-07-16). Every "hybrid vs pure-MPI"
> number produced before 2026-07-16 actually compared hybrid against hybrid-T1.
>
> The sweep therefore builds a **second binary** for the baseline, and it must be
> `-DDENDRO_HYBRID_OMP=OFF -DDENDRO_UNZIP_SPEEDUP=ON`, **not** bare `OFF`. Per
> `dendrolib/CMakeLists.txt:302`, `DENDRO_HYBRID_OMP=ON` *also* silently enables the
> integer-index scatter fast path and the SIMD tensor kernels — both MPI-safe and
> unrelated to threading. A bare-`OFF` baseline gives those up too (measured: unzip
> 193 ms vs 117 ms, **39% slower**) and would **flatter** hybrid by crediting it with
> optimizations it didn't earn. `UNZIP_SPEEDUP=ON` isolates the threading variable.

### The ghost exchange is mostly `Waitall`, and `Waitall` is imbalance

Measured 2026-07-16 (8 cores, depth 9, B/rank 36): **`ghost_wait_s` is 85% of the
exchange at `T=1`**; `pack` + `unpack` together are 15%. So the exchange is dominated
by `MPI_Waitall`, and a rank sits in `Waitall` mainly because its **neighbours are
still computing** — that is load imbalance billed to comm, not wire time. Read a
rising `ghost_pack_wait_s` as an imbalance signal first.

This also means arguments about ghost **volume** are arguing about a minority of the
cost. Two further corrections to the volume intuition, both worth knowing:

- Raising `T` genuinely cuts the job's **aggregate** ghost volume — total ghost bytes
  go as `R·(V/R)^(2/3) ∝ R^(1/3)`, which falls as ranks fall. That intuition is
  correct, but the timers are **per-rank**, and per-rank ghost volume goes as
  `(V/R)^(2/3)`, which **rises ~2.5× at T=4**. Fewer ranks each move more bytes. The
  aggregate win is real and simply invisible to this instrument.
- `pack`/`unpack` are threaded over **neighbour ranks**, not nodes
  (`mesh.tcc:934`, `:1027`), so their parallelism is capped by neighbour count and a
  single large face message cannot be split across threads. Real, but small.

## Reading the CSV

| column | meaning |
|--------|---------|
| `rk_step_s` | **headline** — mean wall time per RK step (lower = better) |
| `rhs_wall_s` `unzip_s` `zip_s` | compute phases (should stay ~flat across the sweep) |
| `ghost_pack_wait_s` | `unzip_wcomm − unzip` = the whole exchange. **Not wire time.** |
| `ghost_pack_s` `ghost_wait_s` `ghost_unpack_s` | that exchange, split. Sum back to `ghost_pack_wait_s` within ~1% **on means** |
| `blocks_per_rank` | **saturation** — keep ≥ 16 or results are noise (script warns if low) |
| `speedup_vs_hybrid_T1` | hybrid-T1 rk_step / this rk_step (>1 → more threads help) |

**Hybrid winning** = flat `rhs_wall_s`, rising speedup.

> ### Timer semantics — read before quoting any phase number
>
> **Quote `rhs_wall_s` for any cross-config comparison.** It comes from `CTX_RHS`,
> which brackets `bssnRHS()` from *outside* the OpenMP region, so it is the whole RHS
> phase as rank wall time.
>
> The JSONL also carries the sub-phase breakdown `deriv_t0`, `rhs_eqn_t0`,
> `rhs_ko_t0`, `bdyc_t0`. These start *inside* the threaded block loop, and
> `profiler_t::start/stop` early-return on worker threads
> (`dendrolib/src/profiler.cpp:28,34`) to avoid racing the shared counter — so each
> accumulates only **thread 0's elapsed time** in that sub-phase. Since thread 0 runs
> concurrently with its peers, that ≈ the rank's wall time for the sub-phase **when
> the threads are balanced**; under imbalance thread 0 need not be the critical path,
> so they understate. `_t0` marks that caveat — it does *not* mean "divide by `T`".
>
> Two traps, both measured (R=2, depth 9, 2026-07-16):
>
> 1. **`t_rhs` excludes `t_deriv`** — it starts *after* the deriv calls
>    (`rhs.cpp:189-205`, then `:221`). Hence `rhs_eqn_t0`, not `rhs_t0`. Reading the
>    old `rhs` key as the whole RHS understates it by ~2×.
> 2. **`t_rhs_ko` is a _subset_ of `t_rhs`** (`profile_params.h:31`), not a sibling.
>    Adding it double-counts by ~12%.
>
> The identity that does hold, and that you should use as a check:
>
> ```
> deriv_t0 + rhs_eqn_t0 + bdyc_t0  ==  rhs_wall        (< 0.2% at T=1 and T=2)
> ```
>
> This is not hypothetical. The old `rhs` key was fed into a
> `rk_step − rhs − uwcomm` subtraction to infer a serial residue; because `rhs`
> was only the equation loop, that residue silently included all of `deriv`, and the
> conclusion drawn from it was wrong.
>
> Profiler barriers (needed to emit the JSONL) inflate the unzip phases, so read
> phase times as **relative**, not absolute.

### Build hygiene: the flag under test must be the ONLY difference

The sweep builds two binaries and pins `DENDRO_USE_NEW_DERIVS=OFF` and
`DVEC_ZERO_ALLOC=ON` on both, and it **refuses to reuse a build tree** whose recorded
`.sweep_flags` don't match what the sweep needs.

That is not paranoia. On 2026-07-16 a hand-built pair silently differed in
`DENDRO_USE_NEW_DERIVS`, which swaps the entire RHS translation unit
(`BSSN_GR/CMakeLists.txt:92` selects `rhs_experimental.cpp` instead of `rhs.cpp`).
The two binaries were running **different physics**, which produced a bogus
"`unzip_scatter_batch` is not bit-exact" verdict and a bogus `all_dg` timing. A
CMake build tree keeps its *own* cached options, so a directory configured in an
earlier session does not take today's defaults — and the binary looks identical from
the outside. If you build variants by hand, diff the caches first:

```bash
diff <(grep -E '^DENDRO_|^BSSN_|^DVEC_' a/CMakeCache.txt|sort) \
     <(grep -E '^DENDRO_|^BSSN_|^DVEC_' b/CMakeCache.txt|sort)
```

And when a bit-exactness gate PASSES, prove it can fail — run it against a
known-different build. Two separate false PASSes were produced that day by diffing
empty files and by diffing a constant awk field.

## Sizing the mesh

The mesh is fixed by the parfile and does **not** grow with node count, so the busiest
config (`T=1`) needs `blocks_per_rank ≥ 16`. The default `q1` BBH grid is ~1100 blocks
— enough for N=1–2, but at N≥4 `T=1` starves (1–2 blocks/rank) and the numbers aren't
reliable. **To saturate N=4/N=8, use the higher-mass-ratio grid** — more AMR blocks with the
same production `ob0` large-block layout that actually shows the hybrid win:

```bash
PARFILE=q8.scaling.par.toml sbatch -A ACCT -p spr --nodes=4 run_hybrid_scaling.sh
```

`q8` is mass ratio 8:1 (`MAXDEPTH=17`); the small BH's deeper refinement gives several× more
blocks. Confirm `blocks_per_rank` on the first run and raise `BSSN_MAXDEPTH` in the parfile
if you need still more.

**Don't** reach for `GRID=uniform` to saturate the hybrid comparison — it's a dead end:
`ob0` uniform makes every block full-size and OOMs, and `ob31` uniform (no fusion) shows **no**
hybrid win (the win is a large-block load-balance effect, absent when all blocks are tiny/equal).
The uniform path stays for other bandwidth experiments; `MEM_CHECK` guards it. `ob31` mechanics:
`blocks ≈ 8^(LEV−3)`, memory-hungry (each tiny block carries a full ghost zone) — `LEV=5` runs
2→8 nodes on one node's RAM, `LEV=7` (~2M blocks) needs ~15+ nodes.

## Roofline mode

`PROFILE=roofline` places the workload against the memory-bandwidth wall (achieved GB/s vs a
STREAM ceiling) — useful for bandwidth-reduction questions like ghost compression. Note the
hybrid *speedup* itself is a large-block load-balance effect, not bandwidth (the `ob31` control
is flat), so this is a bandwidth diagnostic, not an explanation of the hybrid win:

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
| `THREADS_LIST` | `1 2 4` | threads/rank; `1` = hybrid-T1, **not** pure MPI (see above). Keep each a divisor of the per-**socket** core count. |
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

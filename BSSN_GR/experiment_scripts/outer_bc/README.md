---
title: Outer boundary experiments, handoff
updated: 2026-11-23
owner: Nate Garey
advisor: David Van Komen
branch: exp/bssn-outer-bc (experiments only; changes are cherry-picked into experimental once proven)
---

# Outer boundary experiments: per-field Sommerfeld with extrapolated ghosts

> [!summary] The question
> Can the q1 BBH run with the **per-field** outer boundary condition and **extrapolated ghost points** at the outer faces, stay stable through merger and ringdown, and keep the far-field constraint violation well below what we run today (mask 5 + tail)?
> Locally, on a single black hole, it cuts the far-field Hamiltonian constraint about **1000×**. Nobody has run it on the real q1 grid or through merger. That is your project.

## 1. Background

### 1.1 What the outer boundary does
Dendro's grid is a cube (±400 M for q1). At the six faces there are no neighbors, so the evolution equations need a boundary condition. The legacy condition is a radiative (Sommerfeld) condition applied to **every** field $f$:

$$
\partial_t f = -v\left(\frac{x^i}{r}\,\partial_i f + \frac{f - f_\infty}{r}\right),
$$

with $f_\infty$ the flat-space value and $v$ a speed (1 for all fields in the legacy code).

The **per-field** condition (`BSSN_OUTER_BC = 1`) follows Hilditch et al., arXiv:1212.2901, Eqs. (26)–(30), the condition GR-Athena++ also uses. Only $K$, $\tilde\Gamma^i$ and $\tilde A_{ij}$ get a radiative condition, each at its own characteristic speed, and $B^i$ follows $\tilde\Gamma^i$. $\alpha$, $\beta^i$, $\chi$ and $\tilde\gamma_{ij}$ keep their ordinary evolution equations at the face.

**Mask 5 + tail** (`BSSN_OUTER_BC = 3`, `LEGACY_MASK = 5`, `TAIL_MASS = 0.988`) is the current safe choice. It is per-field except that $\alpha$ and $\chi$ stay on the legacy condition, relaxed toward the static puncture profile

$$
f_s = \left(1 + \frac{M}{2r}\right)^{-n},\qquad n = 2 \text{ for } \alpha,\; n = 4 \text{ for } \chi,
\qquad
\partial_t f = -\partial_r (f - f_s) - \frac{f - f_s}{r}.
$$

### 1.2 Why the boundary needs ghosts
Two operators run up to the face:
- **The compact derivatives (E6, A6).** Near a face they switch to one-sided "closure" rows.
- **Kreiss–Oliger (KO) dissipation.** This adds $\sigma\,h^{2r-1}(-1)^{r}\,(D_+D_-)^{r} u$ to each right-hand side. Near the face it also needs special "wall" rows, which are only $O(\sigma h^2)$ accurate.

Earlier tests found two things:
- **The KO wall rows cause most of the far-field error.** Under per-field, they produce almost all of the far-field Hamiltonian violation. Skipping KO on the 4 face rows cut it by 140–1170× locally and 9–53× on BYU.
- **Skipping KO is not a fix.** With no damping at the wall, the **A6** run grew exponentially along the box edges (e-folding time about 3.5 M) and became worse than mask 5 + tail at about 28 M.

**Ghost extrapolation** (`BSSN_OUTER_GHOST_EXTRAP = p`) fills the padding layers outside each physical face with a degree-$p$ Lagrange extrapolation through the face node and the $p$ nodes inside it. At distance $g$ outside the face:

$$
u_{-g} = \sum_{n=0}^{p} w_{g,n}\, u_n,
\qquad
w_{g,n} = \prod_{m \ne n} \frac{-g - m}{n - m}.
$$

Axes are filled in the order x, y, z, so edge and corner ghosts come out as products of the 1D extrapolations. The derivatives then use their **interior** operator everywhere, while the boundary condition still acts on the face nodes.
- With `BSSN_OUTER_GHOST_KO = true`, KO also runs through the ghosts: the face gets real, consistent damping instead of the wall rows.
- With `false`, KO keeps its wall rows.

A desk eigenvalue model (2D/3D boxes with corners) predicts that degree 5 removes A6's corner growth (spectral abscissa 0.33 → 2.5e-8).

> [!important] R065: talk to David before you start
> The A6 arms use the BYU second-derivative operator `BYU_A6_2ND_R065_OP1` ("R065"). Our desk models suggest its wall closure is part of why A6 needs damping at the wall: in a 1D model it has complex wall eigenvalues, and we have no summation-by-parts-type stability proof for it. David will go over what we know about R065 with you. If an A6 result surprises you, R065 is a likely suspect.

## 2. What we already know

**BYU q1, depth 15, ±400 M.** Far-field $C_H$ rms ($r > 300$), as a ratio to mask 5 + tail:

| t / M | current BC | per-field | per-field, KO skipped on 4 face rows |
|---|---|---|---|
| 5 | E6 18× · A6 9.7× | E6 1.5× · A6 2.0× | E6 0.031× · A6 0.019× |
| 20 | E6 9.4× · A6 3.0× | E6 4.9× · A6 7.5× | E6 0.091× · A6 0.11× |
| 30 | E6 6.9× · A6 2.4× | E6 8.1× · A6 12× | E6 0.13× · **A6 1.6× (growing)** |

The reference data is in `ref/`, from jobs 15318687 (current BC), 15335001 (per-field), 15475703 + 15480146 (mask 5 + tail, to 60 M) and 15480423 (KO skip).

**Local single black hole** (M = 1, ±32, depth 6, A6 first derivatives, E6 second derivatives). Far-field $C_H$ rms:

| arm | 6 M | 12 M | 18 M | 24 M |
|---|---|---|---|---|
| per-field | 3.3e-6 | 4.2e-6 | 4.6e-6 | 4.9e-6 |
| + ghosts p5, KO wall rows | 3.2e-6 | 4.1e-6 | 4.5e-6 | 4.9e-6 |
| + ghosts p3, KO wall rows | 3.1e-6 | 4.3e-6 | 4.9e-6 | 5.4e-6 |
| + KO skipped on 4 face rows | 2.4e-9 | 3.6e-9 | 4.4e-9 | 1.3e-8 |
| **+ ghosts p5, KO through ghosts** | **4.1e-9** | **3.7e-9** | **3.0e-9** | **1.0e-8** |
| mask 5 + tail | 3.6e-5 | 4.3e-5 | 5.2e-5 | 5.8e-5 |
| mask 5 + tail, + ghosts p5, KO wall rows | 3.6e-5 | 4.3e-5 | 5.1e-5 | 5.7e-5 |

- **The gain comes from KO running through the ghosts.** Ghosts alone change nothing while KO keeps its wall rows.
- **The late rise is not the BC.** The rise after 24 M shows up in every arm: it is an interior momentum-constraint blow-up near 28 M in this small box.
- **This box can't test the A6 growth.** The KO-skip run also looked fine here, then grew on BYU. **Only the BYU grid can answer the question.**

## 3. Getting set up

1. Code:
   - Dendro-GR at `github.com/dfvankomen/Dendro-GR`, branch `exp/bssn-outer-bc`.
   - dendrolib at `github.com/dfvankomen/Dendro-5.01`, branch `experimental`. It has the BYU A6 operators.
2. Unit test of the ghost fill. Run from the Dendro-GR root; it must print `ALL PASS`:
   ```
   g++ -O2 -std=c++17 -I BSSN_GR/include BSSN_GR/experiment_scripts/outer_bc/test_ghost_extrap.cpp BSSN_GR/src/outer_bc.cpp -o test_ghost_extrap && ./test_ghost_extrap
   ```
3. Initial data: the TwoPunctures file `tp_q001_tpid_sol.bin` for the q1 setup. Ask David for a copy, or generate it with the q1 parfile.
4. Edit the four paths at the top of `run_arms.sbatch` for your account (Eric's allocation). The script builds the solver on first use.
5. Check every run's startup line. A ghost arm must print
   `Outer BC: per-field Sommerfeld (...), ghosts extrapolated at degree 5 (KO through ghosts)`.
   If the ghost part is missing, the parameter was not read. A common cause: a key placed after the `[AEH_PARAMS]` table in the parfile belongs to that table. Keep new keys above it.

## 4. The experiments

Run them in order, and stop to talk with David after each step. One BYU job at a time.

**Step 1: A6 and E6 to 30 M** (`sbatch run_arms.sbatch`, 4 nodes, about 5 h).
- Compare against the references:
  `python3 compare.py mask5tail=ref/A6_mask5tail.dat.gz koskip4=ref/A6_koskip4.dat.gz ghost=<run>/ghost_A6/dat/dgr_ConstraintNorms.dat --png A6.png`, and the same for E6.
- *Pass:* ghost far $C_H$ rms is below mask 5 + tail at every time up to 30 M, **and** the A6 arm does not grow. Over 20–30 M, the slope of $\log C_H$ must be no steeper than mask 5 + tail's. Check `C_HAM_maxabs_far` too, since the KO-skip growth sat on the box edges.
- *Optional control:* add `perfield_A6` (6 nodes). It shows how far your build is from the 15335001 reference. The match won't be bit-for-bit, because the dendrolib snapshot differs.

**Step 2: to 60 M** (`FROM=<step-1 job> TEND=60.0 sbatch restore_arms.sbatch`).
- *Pass:* still below mask 5 + tail, with no growth. The mask 5 + tail reference goes to 60 M.

**Step 3: through merger and ringdown, to 1000 M (production length).**
The kit parfiles already use the production grid: they are `BSSN_GR/pars/q1.par.toml` (depth 15, ±400 M) plus the outer-BC keys. For a merger run you only change the length, switch the horizon finder back on, and write VTUs less often.

- *Arms:* `ghost_E6` and `mask5tail_E6`. E6 first, because it is the production operator and your November E6/E6 run to 1000 M (current BC) is a free third reference. A6 comes after steps 1–2 show its edge growth is gone.
- *Parfiles:* make long copies, so the 30 M ones stay as they are:
  ```
  cd pars
  for a in ghost_E6 mask5tail_E6; do
    sed -e 's/^BSSN_RK_TIME_END = .*/BSSN_RK_TIME_END = 1000.0/' \
        -e 's/^AEH_SOLVER_FREQ = .*/AEH_SOLVER_FREQ = 50/' \
        -e 's/^BSSN_IO_OUTPUT_FREQ = .*/BSSN_IO_OUTPUT_FREQ = 4000/' $a.par.toml > ${a}_long.par.toml
  done
  diff ../../../pars/q1.par.toml ghost_E6_long.par.toml
  ```
  The diff should show only the derivative types, `BSSN_RESTORE_SOLVER`, `BSSN_ASYNC_COMM_K`, the end time, the VTU frequency and the outer-BC keys. VTUs every 4000 steps is about every 16 M. At the kit's 400 an arm would write a few hundred GB by 1000 M (a rough guess from 9.4 GB per 30 M).
- *First piece* (4 nodes per arm, as in production):
  ```
  ARMS="ghost_E6_long mask5tail_E6_long" NPA=4 sbatch --nodes=8 run_arms.sbatch
  ```
  Five hours at ~22 M/h gets to about 100 M, so this piece also covers steps 1 and 2 for E6. Check it with `compare.py` against `ref/E6_*` at 30 and 60 M before going on, and talk to David.
- *The rest:* continue from the last job's checkpoints, one piece at a time:
  ```
  FROM=<last job id> TEND=1000.0 ARMS="ghost_E6_long mask5tail_E6_long" NPA=4 sbatch --nodes=8 restore_arms.sbatch
  ```
  `FROM` is always the job just before, not the first one. Expect about 7 pieces (roughly 18 h to merger at 22 M/h, then about twice as fast; a guess for the ghost arm). After each piece, check `exit.txt` and the restart line in the `.out`, and copy `dat/` and `bah/` somewhere safe: `nobackup/autodelete` purges old files.
- *What to look at, at these times:* the gauge pulse reaches the face centers at ≈ 280 M, the edges at ≈ 400 M and the corners at ≈ 490 M; junk radiation reaches the faces at ≈ 400 M; merger is at ≈ 390 M and its signal reaches the faces at ≈ 795 M.
  1. **Constraints**: `compare.py` with both arms and your November E6/E6 run, `C_HAM_rms_far` and `C_HAM_maxabs_far`. Does either arm jump at those times?
  2. **Horizons**: merger time and final mass should match between the arms. If they don't, the boundary is reaching the binary. AH1/AH2 failing after merger is expected.
  3. **$\Psi_4$ (2,2)** at r = 100: amplitude and phase, ghost vs mask 5 + tail. Compare phase only for t < 560 M (later it is noise). Every November run shows a ~10× amplitude jump at ~785 M, maybe junk reflecting off the wall. Does it shrink with the new BC?
  4. **The faces**: pull VTU frames near 250, 290, 400, 500, 600 and 800 M with `scripts/select_vtu.py` and look at C_HAM on the outer faces and edges. The mask 5 + tail arm holds α to a one-puncture profile; if it is going to break, it should show up when the gauge pulse arrives.
- *Pass:* the ghost arm stays below mask 5 + tail, with no growth at any of those times, and its $\Psi_4$ differs from mask 5 + tail by less than the resolution error from step 4.
- *Untested:* ghosts next to a face during big remeshes (there is a ~4× mesh cliff at ~588 M). If the ghost arm blows up there, say so first, before anything else.

**Step 4: resolution.** Repeat the step-3 arms one depth higher or lower. This gives a resolution error to compare the boundary effect against.

## 5. What to watch for
- **Noise amplification.** The extrapolation weights grow fast. For $p = 5$, $\sum_n |w_{3,n}| = 1023$ at the third ghost, and corners multiply three of these. Noise at the face can be amplified into the ghosts. Look at the outermost rows in the VTU output if anything grows.
- **Untested at the face:** remeshing, refinement changes next to a face, and the GPU build. The GPU build ignores all of these options.
- **Constraint diagnostics** still use the one-sided derivatives at the face. They are a measurement, not part of the evolution.
- **Blocks that are too small** for the degree (fewer than $p + 1$ nodes across) fall back to the one-sided operators without a message. That never happened at element order 6 with $p \le 6$.

## 6. Files here
| file | what |
|---|---|
| `pars/{ghost,perfield,mask5tail}_{A6,E6}.par.toml` | q1 parfiles that differ only in their outer-BC keys (end time 30 M) |
| `ref/*.dat.gz` | BYU `dgr_ConstraintNorms.dat` from the reference jobs above |
| `run_arms.sbatch` | build if needed, run arms from fresh initial data, 2 nodes per arm |
| `restore_arms.sbatch` | continue arms from their checkpoints; cancels itself if the restart is not exact |
| `compare.py` | constraint norms of several runs at fixed times, with ratios and an optional plot |
| `test_ghost_extrap.cpp` | unit test of the ghost fill |

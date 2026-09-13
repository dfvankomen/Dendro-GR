# Experimental on-mesh Teukolsky initial data (type 14)

Use `BSSN_ID_TYPE = 14` with the CPU `bssnSolver` / uniform time stepper
(`TS_MODE=1`). Types 9 and 13 retain their existing initialization behavior.
The initial octree and its initial refinement iterations sample type 9's analytic
conformal metric. After the final mesh is populated, type 14 computes the
connection and seed Ricci, solves on that mesh, then starts evolution.
Restored checkpoints bypass this initial solve.

```toml
BSSN_ID_TYPE = 14
BSSN_DERIVTYPE_FIRST = "E6"
BSSN_DERIVTYPE_SECOND = "E6"
TEUK_HAM_TOL = 1e-10
TEUK_HAM_MAX_ITER = 3000
TEUK_HAM_VERBOSE = true
```

`BSSN_COMPUTE_CONSTRAINTS` and sixth-order derivatives are required. The example
`uniform.toml` uses element order 6, a uniform level-3 mesh (49³ nodes), domain
`[-10,10]³`, amplitude 0.01, pulse center 5, width 2, and zero evolution time.
A normal type-14 run requires no reference file or runtime Python.

## Infrastructure inspection before implementation

The search covered the application sources **and the fetched Dendro source** in
`build/_deps/dendrolib-src`, including ignored dependency files. Terms included
BiCGStab, CG, GMRES, KSP, PETSc, elliptic, Poisson, Jacobi, multigrid, Newton,
matrix-free and dsolve.

* `LinAlg/include/cg.h` implements distributed CG, but calls
  `fem::operators::poisson::matvec` directly, assumes SPD, and writes VTK output
  during iteration. It has no callback for this nonsymmetric finite-difference
  operator. `FEM/include/operators.h` implements FEM Poisson, not the production
  sixth-order BSSN discretization.
* `FEM/include/feMatrix.h`, `include/oda.h`, and the root CMake configuration
  expose optional PETSc support. With `BUILD_WITH_PETSC=ON`, type 14 uses a PETSc
  `MatShell`, MPI vectors containing owned Dendro nodes, and restarted GMRES
  (restart 100, no preconditioner). Matvecs retain Dendro ghost exchange and
  zip/unzip. No global matrix is assembled.
* PETSc is not installed in the validation environment, whose build has
  `BUILD_WITH_PETSC=OFF`. That configuration uses the pre-existing experimental
  BiCGStab, repaired for collective logging, initial zero residual, final true
  residual verification, and collective positivity checks. The PETSc branch
  has not been compiled or run here.
* `BSSN_GR/src/TPNewton.cpp` contains a private BiCGStab and Newton solver tied
  to TwoPunctures' spectral grid. It is not a Dendro nodal solver.
* NLSigma is a nonlinear sigma-model **time-evolution** application.
  `NLSigma/src/nlsmInvCtx.cpp` launches its forward time stepper; its adjoint
  method is a stub. No reusable distributed Newton/GMRES service was found.
  Other Newton routines concern apparent horizons or local fluid inversion.

## Discretization and conventions

The free metric is the determinant-one `U_SYMGT` returned by `NLTeukData`, with
seed `chi=1`. `CodeGen/teuk_connection.py` generates the contracted Christoffel
expression using the same `dendro.py` inverse metric and second Christoffel
symbols as production evolution. Its checked-in include is evaluated with
production metric derivatives on unzipped blocks and zipped onto owned nodes.
This constructs `Gt` before computing seed Ricci and again on the final data.
The independent diagnostic evaluates `Gt - gt^{jk} Gamma^i_jk` without changing
fields; L2/max refer to its Euclidean vector magnitude.

The ordinary `BSSNCtx::compute_constraint_variables()` / `physical_constraints()`
path provides `bar_R = C_HAM` for `chi=1`, `K=At=0`, and consistent `Gt`. No
separate Ricci implementation is introduced. The matrix-free operator uses

```
L u = gt^{ij} D_i D_j u - Gt^i D_i u - bar_R u/8
u = psi - 1,  L u = bar_R/8,  u_initial = 0
```

It uses the same derivative dispatch (including `BSSN_DERIVS` in new-derivative
builds) and successive first derivatives for mixed terms as production Ricci.
Only nodes on physical outer faces have Dirichlet rows `u=0`. Block and MPI
interfaces use Dendro ghost/intergrid operations.

Tolerance is relative Euclidean residual against the initial RHS. Success
requires a **fresh matvec** to satisfy that tolerance, plus finite positive psi;
nonconvergence aborts before evolution and leaves chi unmodified. On success,
`chi=psi^-4` and `U_SYMGT` is retained. This is exactly `gamma=psi^4 bar_gamma`,
`chi=det(gamma)^(-1/3)`, `gt=chi gamma` for determinant-one free data. Lapse stays
1 as in type 9; `K`, `At`, shift and `B` stay zero.

All printed L2 values are RMS over unique owned nodes, including the physical
boundary, rather than volume-weighted AMR norms or the usual BH-masked output.
The elliptic boundary residual is the Dirichlet-row residual. The physical
Hamiltonian is also evaluated at boundary nodes, where no elliptic PDE row is
imposed. Consequently its full-domain norm need not equal the elliptic residual.
Finite-difference derivatives also do not obey an exact continuum chain rule
under `chi=psi^-4`; their truncation errors remain after a converged linear solve.

## Reproduce the comparison

Build as usual, then from `Dendro-GR`:

```sh
cmake --build build --target bssnSolver -j 4
mpirun -np 2 build/BSSN_GR/bssnSolver validation/teuk14/uniform.toml 1 --teuk-compare
python3 validation/teuk14/run_validation.py
```

`--teuk-compare` is optional and requires fresh type-14 data. It samples the
unchanged type-9 and type-13 initializers into scratch arrays on **exactly the
same final Dendro mesh**, evaluates production constraints, and reports pairwise
conformal/physical Frobenius metric differences and chi differences. It does
not overwrite solved evolution fields. `TEUK_SOLVED_ID_FILE` selects the
reference. In this example the binary is `../teukolsky_solved_id.bin` relative
to the process working directory. The regression runner accepts `--reference`,
`--solver`, `--output`, and `--ranks` and resolves paths explicitly.

The regression checks 49³ and 97³ meshes, zero amplitude, inactive ranks, and a
forced nonconvergence. Its common interior region is the central 75% of each
axis (`[-7.5,7.5]³` here), identical in physical extent for both resolutions.
The default output directory is `/tmp/dendro-teuk14-validation`.

## Validation results

See `results.txt` for the measured diagnostics. The reference binary is the
existing Python 49³ solution in the workspace, sampled at coincident nodes for
49³ and trilinearly interpolated for 97³. All rows below use the same mesh for
the three initial-data types.

| Mesh | Type | Full-domain C_HAM RMS | Interior C_HAM RMS | Gamma RMS |
|---|---|---:|---:|---:|
| 49³ | 9 | 7.72945e-4 | 1.17799e-3 | 5.36052e-4 |
| 49³ | 13 | 1.49522e-3 | 8.21545e-8 | 2.75910e-17 |
| 49³ | 14 | 1.97798e-7 | 8.21564e-8 | 0 |
| 97³ | 9 | 8.58499e-4 | 1.31496e-3 | 1.16164e-4 |
| 97³ | 13 | 1.26184e-3 | 1.66958e-3 | 1.75410e-3 |
| 97³ | 14 | 3.37572e-9 | 5.14159e-9 | 0 |

All three momentum RMS values are zero for all rows. Type 14 converged in
126 / 277 BiCGStab iterations on two ranks, with residual RMS
9.43e-15 / 6.08e-15. Its full-domain Hamiltonian improves by 58.6× when spacing
halves (observed order about 5.87). One-rank and two-rank 49³ results agree to
roundoff; iteration counts can differ slightly with MPI reduction ordering.
The flat case has Hamiltonian RMS 2.04e-14 and psi within 3.4e-15 of 1. A CPU
one-step evolution smoke test also completed and wrote a checkpoint; unrelated
wave-extraction `key ... not found` messages occur with this small test domain.

The matching-node interior check agrees with the Python solution's truncation
scale, while the 97³ comparison exposes the file interpolation error. There is
also a distinct boundary difference: the Python solver fixes an outer
**three-cell band** to psi=1 (`solve_hamiltonian.py`, `active[3:-3]`), whereas
type 14 fixes only the physical boundary. Thus the two psi/chi solutions are
not identical, and type 13's full-domain 49³ error is not interpolation error.
Its chi RMS difference from type 14 is 1.68e-4 on 49³. Matching boundary closures
would be necessary for strict pointwise Python/Dendro solution agreement.

This is the CPU uniform-mesh milestone. The implementation uses normal AMR
zip/unzip and physical-face flags, but full AMR accuracy, refinement-interface
convergence, GPU execution, and PETSc execution are not validated here. No
multigrid or preconditioner is supplied; large refined solves may need one.

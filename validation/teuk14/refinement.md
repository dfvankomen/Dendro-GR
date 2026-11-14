# Reducing unnecessary Teukolsky refinement

## What caused it

The original parameters set both BH refinement radii to zero. The old near-BH
test still accepted a node exactly at the origin (`distance <= radius`) and
forced refinement toward the BH maximum level. With BH levels set to 15, this
created a level-13 central cluster. The checkpoint contained only 64 level-13
octants, but they determined the global timestep. Most octants instead resolved
the wave shell, especially at level 8.

Nonpositive radii now disable near-BH forcing; when both are disabled the BH
overrides are bypassed completely. The global minimum depth still applies.
Positive-radius forcing is retained and regression-tested.

Type-14 initialization now resamples the analytic unsolved seed after each trial
remesh. The next wavelet decision therefore assesses the analytic data that will
be populated on that mesh, instead of interpolation artifacts from earlier
meshes. The Hamiltonian solve still runs only after initialization refinement.
This improves consistency but does not by itself stop mesh-size oscillations.

Function 7 remains unchanged: it uses the minimum tolerance inside
`TEUK.REFINEMENT_R0` and increases by a factor of ten per coordinate length unit
until the maximum tolerance. In this case the strict core is radius 2 and the
maximum tolerance is reached at radius 4, while the initial pulse is at radius
20. A different radial function alone cannot remove the central BH override or
put field indicators with different dimensions on a common error budget.

## Optional cell-scale geometric wavelets

Select the new criterion with:

```toml
BSSN_REFINEMENT_MODE = 6
```

The existing tolerance function, refinement variables, coarsening factor,
minimum depth, and active BH overrides are retained. For each element the
wavelet/Sobolev indicator is converted into a dimensionless error:

| Field | Indicator used |
|---|---|
| lapse, chi, shift, conformal metric | original indicator |
| K, conformal connection Gt, shift driver B, conformal At | h × original indicator |

Here `h=max(dx,dy,dz)` is the local nodal spacing. K, Gt, B and At have
inverse-length dimensions in the fixed Cartesian coordinates used here; h times
their error measures their geometric effect across one cell. Refinement still
follows the computed data as the pulse moves or collapses; no prescribed moving
shell, maximum-level cap, or tolerance relaxation is imposed.

This is a **different error criterion**, not a guarantee of pointwise curvature
accuracy at the original absolute tolerance. It is opt-in. Mode 0 keeps its
original field-error criterion, with the zero-radius bug fixed. Mode 6 must be
checked for convergence for each intended evolution, especially near collapse.

## Measured tests

Four-node, 512-rank runs used the original physical parameters and tolerances;
only the endpoint was shortened to t=0.2 for the refinement trials.

| Initialization | Nodes | Maximum level | Initial dt | Hamiltonian iterations | Solved Hamiltonian RMS |
|---|---:|---:|---:|---:|---:|
| Earlier solver-only fix | 17,124,817 | 13 | 0.00101725 | 5200 | 5.20829e-8 |
| Zero-radius fix | 17,048,965 | 9 | 0.0162760 | 4400 | 5.21984e-8 |
| Plus analytic resampling, modes 0 and 6 | 16,642,237 | 9 | 0.0162760 | 4000 | 6.75411e-8 |

The zero-radius fix removed a timestep bottleneck with essentially unchanged
initial Hamiltonian accuracy. Resampling reduced the initial node count modestly
and changed the discrete initial mesh; it did not remove subsequent growth.

In the short fresh evolutions, both resampled runs reached t=0.162760 at step 10:

| Criterion | Nodes after remesh | Maximum level | Next dt | Evolution wall time through t≈0.2 |
|---|---:|---:|---:|---:|
| Mode 0 | 50,888,785 | 10 | 0.00813802 | 41.36 s |
| Mode 6 | 10,304,725 | 8 | 0.0325521 | 20.65 s |

That is 79.8% fewer nodes at the comparison time. The two runs have different
adaptive timesteps after this remesh; their final times/step counts need not
coincide exactly. These timings exclude initialization/solve/output startup and
are not general performance guarantees.

A second comparison restarted **both methods from the same solved checkpoint**
and enabled constraint output every step. The constraints agreed through step
9, before the remeshing comparison. At step 10, t=0.162760:

| Criterion | Hamiltonian RMS | Momentum x/y RMS | Momentum z RMS |
|---|---:|---:|---:|
| Mode 0 | 1.03594e-3 | 8.01534e-4 | 7.38575e-4 |
| Mode 6 | 6.18312e-4 | 3.34859e-4 | 3.73749e-4 |

This supports the reduction in this short test: the smaller mesh did not inflate
the measured constraints. The mode-6 restart completed successfully. The mode-0
restart reached its endpoint, then rank 7 segfaulted in `ot::dealloc_mpi_ctx`
during shutdown; this separate MPI cleanup issue remains unresolved. The
sampled constraint comparison precedes that failure.

This does **not** validate evolution to t=75. Both methods develop substantially larger constraints after transfers than at solved
initialization; long-time AMR/transfer accuracy remains an open issue.

## Reproduce

Use `adaptive.geometric.toml` for a fresh short mode-6 run in a new working
directory with `cp`, `prof`, and `vtu` subdirectories, inside a suitable MPI
allocation. It preserves the original physics, tolerance, and initialization
passes, sets restore=0 and t_end=0.2, and writes constraints every step.
Change only `BSSN_REFINEMENT_MODE` to 0 for the control. Keep the original
`adaptive.toml` unchanged for the earlier full-input reproducer.

Small MPI regression coverage is in `run_solver_validation.py`: zero radii,
positive radii, geometric-mode flat data, uniform/flat solves, inactive ranks,
and convergence guards. Radius tests set the Sobolev weights to zero, matching
the large input and isolating the radius semantics. Individual cases can be run
with `--cases zero_radius positive_radius geometric_flat`.

With type 14 and TEUK_HAM_VERBOSE enabled, AMR logs the number of split candidates
from wavelets, the minimum-depth condition, and the final flags. These counts
help distinguish data-driven refinement from overrides.

The original `function2Octree` bootstrap still uses the global scalar tolerance
before the spatial tolerance function is applied by initial AMR. The initial
construction mesh in these large tests remains about 173 million nodes. The
reported reductions concern the final initialization and evolution meshes, not
peak bootstrap memory. Replacing that bootstrap error budget requires separate
library/API work and validation.


## Full-run configuration

[Teuk.geometric.toml](../../Teuk.geometric.toml) requests t=75 with mode 6,
fresh initialization, the original physical parameters and tolerances, and
separate `teuk_geometric` output prefixes. It restores the original extraction
frequency of 25 for lower output cost. This is a proposed full-run configuration,
not evidence that the t=75 evolution has been validated. The short input above
is retained for comparison and diagnostics.

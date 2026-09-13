#ifndef DENDRO_GR_ON_MESH_TEUKOLSKY_ID_H
#define DENDRO_GR_ON_MESH_TEUKOLSKY_ID_H

#include <utility>

#include "mesh.h"

namespace bssn {

// Reuses the production symbolic contracted Christoffel definition. With
// initialize=false, returns RMS and max of the Euclidean Gamma constraint.
std::pair<double, double> teukolskyConnectionOnMesh(ot::Mesh&, double**,
                                                    bool initialize);

void compareTeukolskyOnMesh(ot::Mesh&, double** solved);

struct TeukHamiltonianSolveResult {
    unsigned int iterations = 0;
    bool converged          = false;
    double residual_l2      = 0.0;
    double residual_max     = 0.0;
    double psi_min = 1.0, psi_max = 1.0;
};

/** Solve for u=psi-1 on the already-converged Dendro mesh and replace chi.
 * vars and seed_ricci are zipped, ghosted nodal arrays.  The conformal metric
 * is retained. Call teukolskyConnectionOnMesh before constructing seed Ricci.
 * Only active mesh ranks participate; tolerance is relative to the initial
 * residual.
 */
TeukHamiltonianSolveResult solveTeukolskyHamiltonianOnMesh(
    ot::Mesh& mesh, double** vars, const double* seed_ricci, double tolerance,
    unsigned int maximum_iterations, bool verbose);

}  // namespace bssn
#endif

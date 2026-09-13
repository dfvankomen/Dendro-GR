// Optional timestep-zero diagnostic. Reference data are sampled into scratch
// arrays; neither the type-9/type-13 initializers nor evolution data are
// changed.
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "OnMeshTeukolskyID.h"
#include "SolvedTeukolskyID.h"
#include "grUtils.h"
#include "physcon.h"

namespace bssn {
namespace {
std::pair<double, double> norm(
    ot::Mesh& mesh, const double* v,
    const std::vector<unsigned char>* mask = nullptr) {
    double local[2] = {0, 0}, sum = 0, maximum = 0;
    unsigned long long count = 0, total = 0;
    for (unsigned n = mesh.getNodeLocalBegin(); n < mesh.getNodeLocalEnd();
         ++n) {
        if (mask && !(*mask)[n]) continue;
        ++count;
        local[0] += v[n] * v[n];
        local[1] = std::max(local[1], std::abs(v[n]));
    }
    const auto comm = mesh.getMPICommunicator();
    MPI_Allreduce(local, &sum, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(local + 1, &maximum, 1, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(&count, &total, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
    return {std::sqrt(sum / total), maximum};
}
void report(ot::Mesh& mesh, int type, double** fields,
            const std::vector<unsigned char>& core) {
    const unsigned nz = mesh.getDegOfFreedom(),
                   nu = mesh.getDegOfFreedomUnZip();
    std::array<std::vector<double>, BSSN_NUM_VARS> uz;
    std::array<const double*, BSSN_NUM_VARS> in;
    std::array<std::vector<double>, BSSN_CONSTRAINT_NUM_VARS> cons;
    std::array<double*, BSSN_CONSTRAINT_NUM_VARS> out;
    for (unsigned s = 0; s < BSSN_NUM_VARS; ++s) {
        uz[s].resize(nu);
        in[s] = uz[s].data();
        mesh.readFromGhostBegin(fields[s], 1);
        mesh.readFromGhostEnd(fields[s], 1);
        mesh.unzip(fields[s], uz[s].data());
    }
    for (unsigned s = 0; s < BSSN_CONSTRAINT_NUM_VARS; ++s) {
        cons[s].resize(nu);
        out[s] = cons[s].data();
    }
    const Point lo = mesh.getDomainMinPt(), hi = mesh.getDomainMaxPt();
    for (const auto& b : mesh.getLocalBlockList()) {
        const unsigned sz[3] = {b.getAllocationSzX(), b.getAllocationSzY(),
                                b.getAllocationSzZ()};
        const double h[3]    = {b.computeDx(lo, hi), b.computeDy(lo, hi),
                                b.computeDz(lo, hi)};
        const double pmin[3] = {
            GRIDX_TO_X(b.getBlockNode().minX()) - BSSN_PADDING_WIDTH * h[0],
            GRIDY_TO_Y(b.getBlockNode().minY()) - BSSN_PADDING_WIDTH * h[1],
            GRIDZ_TO_Z(b.getBlockNode().minZ()) - BSSN_PADDING_WIDTH * h[2]};
        double pmax[3];
        for (unsigned a = 0; a < 3; ++a) pmax[a] = pmin[a] + (sz[a] - 1) * h[a];
        physical_constraints(out.data(), in.data(), b.getOffset(), pmin, pmax,
                             sz, b.getBlkNodeFlag());
    }
    std::vector<double> zipped(nz);
    const char* names[] = {"C_HAM", "C_MOM0", "C_MOM1", "C_MOM2"};
    for (unsigned s = 0; s < 4; ++s) {
        mesh.zip(out[s], zipped.data());
        const auto n        = norm(mesh, zipped.data());
        const auto interior = norm(mesh, zipped.data(), &core);
        if (!mesh.getMPIRank())
            std::cout << "COMPARE core type " << type << " " << names[s]
                      << " L2/max " << interior.first << " " << interior.second
                      << '\n';
        if (!mesh.getMPIRank())
            std::cout << "COMPARE type " << type << " " << names[s]
                      << " L2/max " << n.first << " " << n.second << '\n';
    }
    const auto gamma = teukolskyConnectionOnMesh(mesh, fields, false);
    double low = 1e300, high = -1e300, glo, ghi;
    for (unsigned n = mesh.getNodeLocalBegin(); n < mesh.getNodeLocalEnd();
         ++n) {
        low  = std::min(low, fields[VAR::U_CHI][n]);
        high = std::max(high, fields[VAR::U_CHI][n]);
    }
    MPI_Allreduce(&low, &glo, 1, MPI_DOUBLE, MPI_MIN,
                  mesh.getMPICommunicator());
    MPI_Allreduce(&high, &ghi, 1, MPI_DOUBLE, MPI_MAX,
                  mesh.getMPICommunicator());
    if (!mesh.getMPIRank())
        std::cout << "COMPARE type " << type << " Gamma L2/max " << gamma.first
                  << " " << gamma.second << "\nCOMPARE type " << type
                  << " chi min/max " << glo << " " << ghi << '\n';
}
}  // namespace
void compareTeukolskyOnMesh(ot::Mesh& mesh, double** solved) {
    if (!mesh.isActive()) return;
    const unsigned nz = mesh.getDegOfFreedom();
    std::vector<unsigned char> core(nz, 0);
    std::array<std::vector<double>, BSSN_NUM_VARS> analytic, file;
    std::array<double*, BSSN_NUM_VARS> a, f;
    for (unsigned s = 0; s < BSSN_NUM_VARS; ++s) {
        analytic[s].resize(nz);
        file[s].resize(nz);
        a[s] = analytic[s].data();
        f[s] = file[s].data();
    }
    const auto& elements = mesh.getAllElements();
    const auto& e2n      = mesh.getE2NMapping();
    const auto& dg       = mesh.getE2NMapping_DG();
    const unsigned order = mesh.getElementOrder(),
                   npe   = mesh.getNumNodesPerElement();
    for (unsigned e = mesh.getElementLocalBegin();
         e < mesh.getElementLocalEnd(); ++e)
        for (unsigned p = 0; p < npe; ++p) {
            const unsigned n = e2n[e * npe + p];
            if (n < mesh.getNodeLocalBegin() || n >= mesh.getNodeLocalEnd())
                continue;
            unsigned owner, i, j, k;
            mesh.dg2eijk(dg[e * npe + p], owner, i, j, k);
            const auto& node = elements[owner];
            const double h =
                double(1u << (m_uiMaxDepth - node.getLevel())) / order;
            const double x = node.getX() + i * h, y = node.getY() + j * h,
                         z    = node.getZ() + k * h;
            const double gmax = double(1u << m_uiMaxDepth);
            core[n] = (x >= gmax / 8 && x <= 7 * gmax / 8 && y >= gmax / 8 &&
                       y <= 7 * gmax / 8 && z >= gmax / 8 && z <= 7 * gmax / 8);
            double av[BSSN_NUM_VARS], fv[BSSN_NUM_VARS];
            NLTeukData(x, y, z, av);
            solvedTeukolskyData(x, y, z, fv);
            for (unsigned s = 0; s < BSSN_NUM_VARS; ++s) {
                a[s][n] = av[s];
                f[s][n] = fv[s];
            }
        }
    if (!mesh.getMPIRank())
        std::cout << std::setprecision(16)
                  << "COMPARE identical owned Dendro nodes; L2 denotes nodal "
                     "RMS; core is central 75% of each domain axis\n";
    report(mesh, 9, a.data(), core);
    report(mesh, 13, f.data(), core);
    report(mesh, 14, solved, core);
    // Frobenius norm includes both off-diagonal entries. Report both conformal
    // and physical metrics, and chi separately, for all three pairs.
    double** sets[]   = {a.data(), f.data(), solved};
    const int types[] = {9, 13, 14};
    std::vector<double> difference(nz);
    for (unsigned x = 0; x < 3; ++x)
        for (unsigned y = x + 1; y < 3; ++y)
            for (unsigned kind = 0; kind < 3; ++kind) {
                for (unsigned n = mesh.getNodeLocalBegin();
                     n < mesh.getNodeLocalEnd(); ++n) {
                    double sum = 0;
                    if (kind == 2)
                        sum = std::pow(
                            sets[x][VAR::U_CHI][n] - sets[y][VAR::U_CHI][n], 2);
                    else
                        for (unsigned s = 0; s < 6; ++s) {
                            double vx = sets[x][VAR::U_SYMGT0 + s][n],
                                   vy = sets[y][VAR::U_SYMGT0 + s][n];
                            if (kind == 1) {
                                vx /= sets[x][VAR::U_CHI][n];
                                vy /= sets[y][VAR::U_CHI][n];
                            }
                            sum += (s == 1 || s == 2 || s == 4 ? 2 : 1) *
                                   (vx - vy) * (vx - vy);
                        }
                    difference[n] = std::sqrt(sum);
                }
                const auto n = norm(mesh, difference.data());
                if (!mesh.getMPIRank())
                    std::cout << "COMPARE difference " << types[x] << "-"
                              << types[y] << " "
                              << (kind == 0   ? "conformal metric"
                                  : kind == 1 ? "physical metric"
                                              : "chi")
                              << " L2/max " << n.first << " " << n.second
                              << '\n';
            }
}
}  // namespace bssn

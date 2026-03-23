#ifndef BSSN_SHELL_INTERP_H
#define BSSN_SHELL_INTERP_H

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "mesh.h"
#include "grDef.h"

namespace bssn {
namespace shell_interp {

inline double smoothstep5(double s) {
    if (s <= 0.0) return 0.0;
    if (s >= 1.0) return 1.0;
    return s * s * s * (10.0 + s * (-15.0 + 6.0 * s));
}

struct ShellSampleInfo {
    bool in_shell;
    double s;
    double xin[3];
    double xout[3];
};

inline ShellSampleInfo spherical_shell_map(double x, double y, double z,
                                           double r0, double width) {
    ShellSampleInfo out{};
    out.in_shell = false;
    out.s = 0.0;

    const double rin  = r0 - 0.5 * width;
    const double rout = r0 + 0.5 * width;

    const double r = std::sqrt(x * x + y * y + z * z);
    if (r <= 1.0e-14) return out;
    if (r < rin || r > rout) return out;

    const double nx = x / r;
    const double ny = y / r;
    const double nz = z / r;

    out.in_shell = true;
    out.s = (r - rin) / (rout - rin);

    out.xin[0] = rin * nx;
    out.xin[1] = rin * ny;
    out.xin[2] = rin * nz;

    out.xout[0] = rout * nx;
    out.xout[1] = rout * ny;
    out.xout[2] = rout * nz;

    return out;
}

inline void get_node_xyz(const ot::Mesh* pMesh, unsigned int node,
                         double& X, double& Y, double& Z) {
    const ot::TreeNode* allElements = pMesh->getAllElements().data();
    const auto& cgToDg              = pMesh->getCG2DGMap();
    const unsigned int eleOrder     = pMesh->getElementOrder();

    unsigned int ownerID, ii_x, jj_y, kk_z;
    const unsigned int lookup = cgToDg[node];

    pMesh->dg2eijk(lookup, ownerID, ii_x, jj_y, kk_z);

    const ot::TreeNode& tmpOct = allElements[ownerID];
    const double hx = (tmpOct.maxX() - tmpOct.minX()) / ((double)eleOrder);

    const double x = tmpOct.minX() + ii_x * hx;
    const double y = tmpOct.minY() + jj_y * hx;
    const double z = tmpOct.minZ() + kk_z * hx;

    X = GRIDX_TO_X(x);
    Y = GRIDY_TO_Y(y);
    Z = GRIDZ_TO_Z(z);
}

struct NodeCache {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<unsigned int> nodes;  // actual mesh node ids
};

inline NodeCache build_node_cache(const ot::Mesh* pMesh) {
    NodeCache cache;

    const unsigned int n0 = pMesh->getNodeLocalBegin();
    const unsigned int n1 = pMesh->getNodeLocalEnd();
    const unsigned int n  = n1 - n0;

    cache.x.resize(n);
    cache.y.resize(n);
    cache.z.resize(n);
    cache.nodes.resize(n);

    unsigned int idx = 0;
    for (unsigned int node = n0; node < n1; node++, idx++) {
        cache.nodes[idx] = node;
        get_node_xyz(pMesh, node, cache.x[idx], cache.y[idx], cache.z[idx]);
    }

    return cache;
}

inline double sample_var_nearest(const ot::Mesh* pMesh,
                                 const DendroScalar* var,
                                 double xs, double ys, double zs) {
    double best_r2  = std::numeric_limits<double>::max();
    double best_val = 0.0;

    for (unsigned int node = pMesh->getNodeLocalBegin();
         node < pMesh->getNodeLocalEnd(); node++) {

        double X, Y, Z;
        get_node_xyz(pMesh, node, X, Y, Z);

        const double dx = X - xs;
        const double dy = Y - ys;
        const double dz = Z - zs;
        const double r2 = dx * dx + dy * dy + dz * dz;

        if (r2 < best_r2) {
            best_r2  = r2;
            best_val = var[node];
        }
    }

    return best_val;
}

inline double sample_var_idw(const NodeCache& cache,
                             const DendroScalar* var,
                             double xs, double ys, double zs,
                             unsigned int k = 8,
                             double power = 2.0) {
    if (cache.nodes.empty()) return 0.0;
    if (k == 0) k = 1;
    if (k > cache.nodes.size()) k = static_cast<unsigned int>(cache.nodes.size());

    std::vector<std::pair<double, unsigned int>> best;
    best.reserve(k);

    for (unsigned int i = 0; i < cache.nodes.size(); i++) {
        const double dx = cache.x[i] - xs;
        const double dy = cache.y[i] - ys;
        const double dz = cache.z[i] - zs;
        const double r2 = dx * dx + dy * dy + dz * dz;

        // exact hit
        if (r2 < 1.0e-30) return var[cache.nodes[i]];

        if (best.size() < k) {
            best.emplace_back(r2, i);
            if (best.size() == k) {
                std::sort(best.begin(), best.end(),
                          [](const auto& a, const auto& b) {
                              return a.first < b.first;
                          });
            }
        } else if (r2 < best.back().first) {
            best.back() = {r2, i};
            std::sort(best.begin(), best.end(),
                      [](const auto& a, const auto& b) {
                          return a.first < b.first;
                      });
        }
    }

    double wsum = 0.0;
    double vsum = 0.0;

    for (const auto& entry : best) {
        const double r2 = entry.first;
        const unsigned int idx = entry.second;
        const double w = 1.0 / std::pow(std::sqrt(r2), power);

        wsum += w;
        vsum += w * var[cache.nodes[idx]];
    }

    return (wsum > 0.0) ? (vsum / wsum) : var[cache.nodes[best[0].second]];
}

}  // namespace shell_interp
}  // namespace bssn

#endif
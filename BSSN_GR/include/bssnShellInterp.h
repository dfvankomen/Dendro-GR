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

struct ShellSampleInfo {
    bool in_shell;
    double s;
    double r;
    double rin;
    double rout;
    double n[3];
    double xin[3];
    double xout[3];
};

inline ShellSampleInfo spherical_shell_map(double x, double y, double z,
                                           double r0, double width) {
    ShellSampleInfo out{};
    out.in_shell = false;
    out.s        = 0.0;
    out.r        = 0.0;
    out.rin      = r0 - 0.5 * width;
    out.rout     = r0 + 0.5 * width;
    out.n[0]     = 0.0;
    out.n[1]     = 0.0;
    out.n[2]     = 0.0;

    const double r = std::sqrt(x * x + y * y + z * z);
    out.r = r;

    if (r <= 1.0e-14) return out;
    if (r < out.rin || r > out.rout) return out;

    const double nx = x / r;
    const double ny = y / r;
    const double nz = z / r;

    out.in_shell = true;
    out.s        = (r - out.rin) / (out.rout - out.rin);

    out.n[0] = nx;
    out.n[1] = ny;
    out.n[2] = nz;

    out.xin[0] = out.rin * nx;
    out.xin[1] = out.rin * ny;
    out.xin[2] = out.rin * nz;

    out.xout[0] = out.rout * nx;
    out.xout[1] = out.rout * ny;
    out.xout[2] = out.rout * nz;

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
    std::vector<unsigned int> nodes;
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

inline double sample_on_ray_idw(const NodeCache& cache,
                                const DendroScalar* var,
                                const double n[3], double rr,
                                unsigned int k = 8,
                                double power = 2.0) {
    return sample_var_idw(cache, var,
                          rr * n[0], rr * n[1], rr * n[2],
                          k, power);
}

inline double deriv_inner_2nd_order(double fm2, double fm1, double f0, double h) {
    return (3.0 * f0 - 4.0 * fm1 + fm2) / (2.0 * h);
}

inline double deriv_outer_2nd_order(double f0, double fp1, double fp2, double h) {
    return (-3.0 * f0 + 4.0 * fp1 - fp2) / (2.0 * h);
}

inline double hermite_cubic(double fin, double dfin,
                            double fout, double dfout,
                            double rin, double rout, double r) {
    const double dr = rout - rin;
    if (dr <= 0.0) return fin;

    double t = (r - rin) / dr;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;

    const double t2 = t * t;
    const double t3 = t2 * t;

    const double h00 =  2.0 * t3 - 3.0 * t2 + 1.0;
    const double h10 =        t3 - 2.0 * t2 + t;
    const double h01 = -2.0 * t3 + 3.0 * t2;
    const double h11 =        t3 -       t2;

    return h00 * fin + h10 * dr * dfin
         + h01 * fout + h11 * dr * dfout;
}

inline double sample_var_cubic_hermite_shell(const NodeCache& cache,
                                             const DendroScalar* var,
                                             const ShellSampleInfo& sh,
                                             double deriv_h,
                                             unsigned int k = 8,
                                             double power = 2.0) {
    const double rin  = sh.rin;
    const double rout = sh.rout;
    const double r    = sh.r;

    const double width = rout - rin;
    if (width <= 0.0) return 0.0;

    double h = deriv_h;
    if (h <= 0.0) h = 0.1 * width;

    const double fin = sample_on_ray_idw(cache, var, sh.n, rin,  k, power);
    const double fout = sample_on_ray_idw(cache, var, sh.n, rout, k, power);

    // Use trusted points outside the shell to estimate boundary derivatives
    double rin_m1 = rin - h;
    double rin_m2 = rin - 2.0 * h;
    double rout_p1 = rout + h;
    double rout_p2 = rout + 2.0 * h;

    // keep radii positive near the origin
    if (rin_m1 <= 1.0e-12) rin_m1 = 0.5 * rin;
    if (rin_m2 <= 1.0e-12) rin_m2 = 0.25 * rin;

    // If we had to clip, use the actual spacing that resulted
    const double hin1 = rin - rin_m1;
    const double hin2 = rin_m1 - rin_m2;

    double dfin = 0.0;
    if (std::abs(hin1 - hin2) / std::max(hin1, 1.0e-14) < 1.0e-8) {
        const double fm1 = sample_on_ray_idw(cache, var, sh.n, rin_m1, k, power);
        const double fm2 = sample_on_ray_idw(cache, var, sh.n, rin_m2, k, power);
        dfin = deriv_inner_2nd_order(fm2, fm1, fin, hin1);
    } else {
        const double fm1 = sample_on_ray_idw(cache, var, sh.n, rin_m1, k, power);
        dfin = (fin - fm1) / std::max(rin - rin_m1, 1.0e-14);
    }

    const double fp1 = sample_on_ray_idw(cache, var, sh.n, rout_p1, k, power);
    const double fp2 = sample_on_ray_idw(cache, var, sh.n, rout_p2, k, power);
    const double dfout = deriv_outer_2nd_order(fout, fp1, fp2, h);

    return hermite_cubic(fin, dfin, fout, dfout, rin, rout, r);
}

}  // namespace shell_interp
}  // namespace bssn

#endif
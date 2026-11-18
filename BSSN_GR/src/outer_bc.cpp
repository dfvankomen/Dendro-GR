// Optional outer boundary conditions (BSSN_OUTER_BC); see outer_bc.h.
#include "outer_bc.h"

#include <cmath>

namespace bssn {

namespace {

struct FaceNode {
    unsigned int pp;  // index into the padded block
    double x[3];
    double r;
    int side[3];  // -1 / +1 on a low / high physical face along that axis
};

// Calls fn(FaceNode) once for every node on a physical face of the block.
// Bit 2a of bflag marks the low face along axis a, bit 2a + 1 the high face.
template <typename Fn>
void for_each_face_node(const double *pmin, const double *pmax,
                        const unsigned int *sz, unsigned int pad,
                        unsigned int bflag, Fn &&fn) {
    if (bflag == 0u) return;
    const unsigned int n[3] = {sz[0], sz[1], sz[2]};
    double h[3];
    bool lo[3], hi[3];
    for (unsigned int a = 0; a < 3; a++) {
        h[a]  = (pmax[a] - pmin[a]) / (double)(n[a] - 1);
        lo[a] = bflag & (1u << (2 * a));
        hi[a] = bflag & (1u << (2 * a + 1));
    }
    for (unsigned int k = pad; k < n[2] - pad; k++)
        for (unsigned int j = pad; j < n[1] - pad; j++) {
            const bool jk_face =
                (lo[1] && j == pad) || (hi[1] && j == n[1] - pad - 1) ||
                (lo[2] && k == pad) || (hi[2] && k == n[2] - pad - 1);
            for (unsigned int i = pad; i < n[0] - pad; i++) {
                const bool i_face =
                    (lo[0] && i == pad) || (hi[0] && i == n[0] - pad - 1);
                if (!jk_face && !i_face) {
                    // interior of the row: jump to the high face, if any
                    if (hi[0] && i < n[0] - pad - 1) {
                        i = n[0] - pad - 2;
                        continue;
                    }
                    break;
                }
                const unsigned int idx[3] = {i, j, k};
                FaceNode p;
                p.pp = i + n[0] * (j + n[1] * k);
                for (unsigned int a = 0; a < 3; a++) {
                    p.x[a]    = pmin[a] + idx[a] * h[a];
                    p.side[a] = (lo[a] && idx[a] == pad)              ? -1
                                : (hi[a] && idx[a] == n[a] - pad - 1) ? +1
                                                                      : 0;
                }
                p.r = std::sqrt(p.x[0] * p.x[0] + p.x[1] * p.x[1] +
                                p.x[2] * p.x[2]);
                fn(p);
            }
        }
}

// d_axis f at pp: the given derivative g, or with deriv2 on a face plane the
// second-order one-sided stencil pointing into the block.
inline double d_axis(const double *f, const double *g, unsigned int pp,
                     int side, unsigned int stride, double h, bool deriv2) {
    if (!deriv2 || side == 0) return g[pp];
    if (side < 0)
        return (-3.0 * f[pp] + 4.0 * f[pp + stride] - f[pp + 2 * stride]) /
               (2.0 * h);
    return (3.0 * f[pp] - 4.0 * f[pp - stride] + f[pp - 2 * stride]) /
           (2.0 * h);
}

// 1+log slicing: sqrt(mu_L) alpha with mu_L = 2/alpha
inline double gauge_speed(double alpha) { return std::sqrt(2.0 * alpha); }

// Gamma-driver: the paper's mu_S alpha^2 is Dendro's
// 3/4 (lambda_f[0] + lambda_f[1] alpha). Tangential speed sqrt(mu_S) alpha,
// normal speed 2/sqrt(3) sqrt(mu_S) alpha.
inline double shift_speed_t(double alpha, const double lambda_f[2]) {
    return std::sqrt(0.75 * (lambda_f[0] + lambda_f[1] * alpha));
}
inline double shift_speed_n(double alpha, const double lambda_f[2]) {
    return 2.0 / std::sqrt(3.0) * shift_speed_t(alpha, lambda_f);
}

}  // namespace

void outer_bc_hilditch(const OuterBCFields &u, const OuterBCGrads &g,
                       const OuterBCRhs &r, const double *pmin,
                       const double *pmax, const unsigned int *sz,
                       unsigned int pad, unsigned int bflag,
                       const double lambda_f[2], bool deriv2) {
    const unsigned int stride[3] = {1u, sz[0], sz[0] * sz[1]};
    double h[3];
    for (unsigned int a = 0; a < 3; a++)
        h[a] = (pmax[a] - pmin[a]) / (double)(sz[a] - 1);

    for_each_face_node(pmin, pmax, sz, pad, bflag, [&](const FaceNode &p) {
        const unsigned int pp = p.pp;
        const double inv_r    = 1.0 / p.r;
        const double s[3] = {p.x[0] * inv_r, p.x[1] * inv_r, p.x[2] * inv_r};
        const double al   = u.alpha[pp];
        const double b[3] = {u.beta[0][pp], u.beta[1][pp], u.beta[2][pp]};

        // ds = s.d f, adv = beta.d f
        auto ds_adv       = [&](const double *f, const double *const gf[3],
                                double &ds, double &adv) {
            ds = adv = 0.0;
            for (unsigned int a = 0; a < 3; a++) {
                const double da =
                    d_axis(f, gf[a], pp, p.side[a], stride[a], h[a], deriv2);
                ds += s[a] * da;
                adv += b[a] * da;
            }
        };
        double ds, adv;

        // Eq. (26)
        ds_adv(u.K, g.K, ds, adv);
        r.K[pp] = -gauge_speed(al) * (ds + u.K[pp] * inv_r) + adv;

        // Eq. (30), without the paper's trace projection
        for (unsigned int c = 0; c < 6; c++) {
            ds_adv(u.At[c], g.At[c], ds, adv);
            r.At[c][pp] = -al * (ds + u.At[c][pp] * inv_r) + adv;
        }

        // Eqs. (27)-(28): normal and tangential parts of Gt
        const double vT = shift_speed_t(al, lambda_f);
        const double vL = shift_speed_n(al, lambda_f);
        double G[3], dsG[3], advG[3];
        for (unsigned int c = 0; c < 3; c++) {
            G[c] = u.Gt[c][pp];
            ds_adv(u.Gt[c], g.Gt[c], dsG[c], advG[c]);
        }
        const double Gs    = s[0] * G[0] + s[1] * G[1] + s[2] * G[2];
        const double dGs   = s[0] * dsG[0] + s[1] * dsG[1] + s[2] * dsG[2];
        const double aGs   = s[0] * advG[0] + s[1] * advG[1] + s[2] * advG[2];
        const double rhs_s = -vL * (dGs + Gs * inv_r) + aGs;
        for (unsigned int c = 0; c < 3; c++) {
            const double GA  = G[c] - s[c] * Gs;
            const double dGA = dsG[c] - s[c] * dGs;
            const double aGA = advG[c] - s[c] * aGs;
            const double gt_bc =
                s[c] * rhs_s + (-vT * (dGA + GA * inv_r) + aGA);
            // The driver is d_t B = d_t Gt + (terms without d_t Gt): swap the
            // interior d_t Gt for the boundary one.
            r.B[c][pp] += gt_bc - r.Gt[c][pp];
            r.Gt[c][pp] = gt_bc;
        }
    });
}

void outer_bc_speed(double *f_rhs, const double *f, const double *const df[3],
                    const double *alpha, const double *pmin, const double *pmax,
                    const unsigned int *sz, unsigned int pad,
                    unsigned int bflag, double falloff, double f_inf,
                    BCSpeed kind, const double lambda_f[2]) {
    for_each_face_node(pmin, pmax, sz, pad, bflag, [&](const FaceNode &p) {
        const unsigned int pp = p.pp;
        const double al       = alpha[pp];
        const double v = kind == BCSpeed::Gauge   ? gauge_speed(al)
                         : kind == BCSpeed::Shift ? shift_speed_n(al, lambda_f)
                                                  : al;
        f_rhs[pp]      = -v * (1.0 / p.r) *
                         (p.x[0] * df[0][pp] + p.x[1] * df[1][pp] +
                          p.x[2] * df[2][pp] + falloff * (f[pp] - f_inf));
    });
}

void outer_bc_tail(double *f_rhs, const double *pmin, const double *pmax,
                   const unsigned int *sz, unsigned int pad, unsigned int bflag,
                   double M, double n) {
    for_each_face_node(pmin, pmax, sz, pad, bflag, [&](const FaceNode &p) {
        const double q   = 1.0 + 0.5 * M / p.r;
        const double fs  = std::pow(q, -n);
        const double dfs = n * 0.5 * M / (p.r * p.r) * std::pow(q, -n - 1.0);
        f_rhs[p.pp] += dfs + (fs - 1.0) / p.r;
    });
}

}  // namespace bssn

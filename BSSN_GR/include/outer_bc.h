// Optional outer boundary conditions for BSSN, selected by BSSN_OUTER_BC
// (rhs_experimental.cpp builds only). The legacy bssn_bcs() is mode 0.
//
// Mode 1 is the per-field Sommerfeld condition of D. Hilditch et al.,
// "Compact binary evolutions with the Z4c formulation", Phys. Rev. D 88,
// 084057 (2013), arXiv:1212.2901, Eqs. (26)-(30), applied to the BSSN
// variables with Theta dropped. GR-Athena++ (arXiv:2101.08289, Secs. 3.2,
// 3.4) uses the same conditions on an octree cube. Only K, Gt and At get a
// radiative condition; alpha, beta, chi and gt keep their interior equations
// at the face, and B follows Gt.
#pragma once

namespace bssn {

struct OuterBCFields {
    const double *alpha, *beta[3], *K, *Gt[3], *At[6];
};
struct OuterBCGrads {  // first derivatives, [field][axis]
    const double *K[3], *Gt[3][3], *At[6][3];
};
struct OuterBCRhs {  // overwritten on face nodes
    double *K, *Gt[3], *At[6], *B[3];
};

// Mode 1: Eqs. (26)-(28), (30) of arXiv:1212.2901 on the physical-face nodes
// of one block. lambda_f is BSSN_LAMBDA_F (shift driver
// d_t beta = 3/4 (lambda_f[0] + lambda_f[1] alpha) B + ...). With deriv2 the
// normal derivative on face nodes uses a second-order one-sided stencil
// (GR-Athena++ instead uses second-order centered differences on
// extrapolated ghosts, which Dendro does not have at physical faces).
void outer_bc_hilditch(const OuterBCFields &u, const OuterBCGrads &g,
                       const OuterBCRhs &r, const double *pmin,
                       const double *pmax, const unsigned int *sz,
                       unsigned int pad, unsigned int bflag,
                       const double lambda_f[2], bool deriv2);

// Mode 2: the legacy condition d_t f = -v [s.d f + falloff (f - f_inf) / r]
// at the field's own characteristic speed v (the speeds of Eqs. (26)-(30);
// Gt, beta and B use the normal shift speed).
enum class BCSpeed { Gauge, Shift, Light };
void outer_bc_speed(double *f_rhs, const double *f, const double *const df[3],
                    const double *alpha, const double *pmin, const double *pmax,
                    const unsigned int *sz, unsigned int pad,
                    unsigned int bflag, double falloff, double f_inf,
                    BCSpeed kind, const double lambda_f[2]);

// BSSN_OUTER_BC_TAIL_MASS: adds d_r f_s + (f_s - 1)/r on face nodes after the
// legacy condition, f_s = (1 + M/2r)^-n, so the condition becomes
// d_t f = -d_r (f - f_s) - (f - f_s)/r and holds the static puncture tail
// (n = 2 for alpha, 4 for chi). The legacy condition alone is static only for
// f = 1 + c/r and drives the M^2/r^2 part of the tail.
void outer_bc_tail(double *f_rhs, const double *pmin, const double *pmax,
                   const unsigned int *sz, unsigned int pad, unsigned int bflag,
                   double M, double n);

// BSSN_OUTER_GHOST_EXTRAP: fills the padding beyond the physical faces of a
// block by degree-p Lagrange extrapolation through the face node and the p
// nodes inside it. Axes go x, y, z, so edge and corner ghosts are tensor-
// product extrapolations. Call only where outer_ghost_extrap_fits().
constexpr unsigned int kMaxGhostDegree = 8, kMaxGhostPad = 8;
bool outer_ghost_extrap_fits(const unsigned int *sz, unsigned int pad,
                             unsigned int p);
void outer_ghost_extrap(double *f, const unsigned int *sz, unsigned int pad,
                        unsigned int bflag, unsigned int p);

}  // namespace bssn

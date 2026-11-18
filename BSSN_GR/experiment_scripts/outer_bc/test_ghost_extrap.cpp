// Standalone check of bssn::outer_ghost_extrap (BSSN_OUTER_GHOST_EXTRAP).
//
// Build and run from the repository root:
//   g++ -O2 -std=c++17 -I BSSN_GR/include
//       BSSN_GR/experiment_scripts/outer_bc/test_ghost_extrap.cpp
//       BSSN_GR/src/outer_bc.cpp -o test_ghost_extrap && ./test_ghost_extrap
//
// Passes if degree-p extrapolation reproduces a degree-p polynomial in every
// ghost (faces, edges, corners) for all 63 face combinations, leaves every
// other node untouched, is not exact one degree higher, and rejects blocks
// with too few nodes.
#include <cmath>
#include <cstdio>
#include <vector>

#include "outer_bc.h"

namespace {

constexpr unsigned int kPad = 3;

struct Block {
    unsigned int sz[3];
    std::vector<double> f;
    explicit Block(unsigned int nx, unsigned int ny, unsigned int nz)
        : sz{nx, ny, nz}, f((size_t)nx * ny * nz) {}
    double &at(unsigned int i, unsigned int j, unsigned int k) {
        return f[i + (size_t)sz[0] * (j + (size_t)sz[1] * k)];
    }
};

// Degree p along each axis, plus one mixed term.
double poly(double x, double y, double z, unsigned int p) {
    double v = 0.3 + 0.05 * x * y * z;
    for (unsigned int a = 1; a <= p; a++)
        v += 0.7 / a * std::pow(x, a) - 0.2 * a * std::pow(y, a) +
             (a % 2 ? 0.11 : -0.11) * std::pow(z, a);
    return v;
}

bool in_physical_ghost(const Block &b, unsigned int bflag,
                       const unsigned int idx[3]) {
    for (unsigned int axis = 0; axis < 3; axis++) {
        if ((bflag & (1u << (2 * axis))) && idx[axis] < kPad) return true;
        if ((bflag & (1u << (2 * axis + 1))) && idx[axis] >= b.sz[axis] - kPad)
            return true;
    }
    return false;
}

// Roundoff bound for a corner ghost: eps * (sum |w| at the outermost ghost)^3.
double roundoff_bound(unsigned int p) {
    double l1 = 0.0;
    for (unsigned int n = 0; n <= p; n++) {
        double w = 1.0;
        for (unsigned int m = 0; m <= p; m++)
            if (m != n) w *= (-(double)kPad - m) / ((double)n - m);
        l1 += std::fabs(w);
    }
    return 8.0 * 2.2e-16 * l1 * l1 * l1;
}

// Max |f - exact| over the block, infinite if any node is NaN.
double max_error(Block &b, unsigned int p) {
    double err = 0.0;
    for (unsigned int k = 0; k < b.sz[2]; k++)
        for (unsigned int j = 0; j < b.sz[1]; j++)
            for (unsigned int i = 0; i < b.sz[0]; i++) {
                const double d = std::fabs(
                    b.at(i, j, k) - poly(0.1 * i, 0.13 * j, 0.09 * k, p));
                if (std::isnan(d)) return INFINITY;
                err = std::fmax(err, d);
            }
    return err;
}

int check_exact_on_polynomials() {
    int fails = 0;
    for (unsigned int p = 1; p <= 6; p++)
        for (unsigned int n : {7u, 9u})
            for (unsigned int bflag = 1; bflag < 64; bflag++) {
                Block b(n + 2 * kPad, n + 1 + 2 * kPad, n + 2 + 2 * kPad);
                for (unsigned int k = 0; k < b.sz[2]; k++)
                    for (unsigned int j = 0; j < b.sz[1]; j++)
                        for (unsigned int i = 0; i < b.sz[0]; i++) {
                            const unsigned int idx[3] = {i, j, k};
                            b.at(i, j, k) =
                                in_physical_ghost(b, bflag, idx)
                                    ? NAN
                                    : poly(0.1 * i, 0.13 * j, 0.09 * k, p);
                        }
                if (!bssn::outer_ghost_extrap_fits(b.sz, kPad, p)) {
                    printf("FAIL: p %u, n %u rejected\n", p, n);
                    fails++;
                    continue;
                }
                bssn::outer_ghost_extrap(b.f.data(), b.sz, kPad, bflag, p);
                const double err = max_error(b, p);
                if (!(err <= roundoff_bound(p) * 10.0)) {
                    printf("FAIL: p %u, n %u, bflag %u, error %.3e\n", p, n,
                           bflag, err);
                    fails++;
                }
            }
    return fails;
}

int check_not_exact_one_degree_higher() {
    Block b(13, 13, 13);
    for (unsigned int k = 0; k < 13; k++)
        for (unsigned int j = 0; j < 13; j++)
            for (unsigned int i = 0; i < 13; i++)
                b.at(i, j, k) = std::pow(0.1 * i, 6);
    bssn::outer_ghost_extrap(b.f.data(), b.sz, kPad, 1u, 5);
    const double err = std::fabs(b.at(0, 6, 6) - std::pow(-0.3, 6));
    if (err > 1e-8) return 0;
    printf("FAIL: degree 5 exact on x^6 (error %.3e)\n", err);
    return 1;
}

int check_small_block_rejected() {
    const unsigned int sz[3] = {2 * kPad + 5, 13, 13};
    if (!bssn::outer_ghost_extrap_fits(sz, kPad, 5)) return 0;
    printf("FAIL: 5 nodes across accepted for degree 5\n");
    return 1;
}

}  // namespace

int main() {
    const int fails = check_exact_on_polynomials() +
                      check_not_exact_one_degree_higher() +
                      check_small_block_rejected();
    printf("%s (%d failures)\n", fails ? "FAILED" : "ALL PASS", fails);
    return fails != 0;
}

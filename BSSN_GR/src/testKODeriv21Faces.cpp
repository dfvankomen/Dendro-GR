// Oracle: the closure as coded. On a 13^3 random block, ko_deriv21_y/z rows
// 2, 3 (low) or n-4, n-3 (high) equal (u0 - 2 u1 + u2) / (4 h) on the nearest
// three points under the block's own face flag (DOWN/UP, BACK/FRONT), and the
// 5-point interior stencil under an x flag (LEFT). Catches: flag mix-ups, not
// closure accuracy. Usage: testKODeriv21Faces

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "derivs.h"

int main() {
    const int n = 13, nx = n, ny = n;  // IDX() reads nx, ny
    const unsigned int sz[3] = {n, n, n};
    const double h           = 0.1;
    std::vector<double> u(n * n * n), Du(n * n * n);
    std::mt19937 gen(3);
    std::uniform_real_distribution<double> r(-1.0, 1.0);
    for (double& v : u) v = r(gen);

    struct Case {
        const char* name;
        unsigned int dir;
        int rows[2];
        int base;  // first point of the closure stencil; -1 for interior
    };
    int bad = 0;
    for (int axis = 1; axis <= 2; axis++) {
        const int stride   = axis == 1 ? nx : nx * ny;
        const Case cases[] = {{axis == 1 ? "DOWN" : "BACK",
                               axis == 1 ? OCT_DIR_DOWN : OCT_DIR_BACK,
                               {2, 3},
                               2},
                              {axis == 1 ? "UP" : "FRONT",
                               axis == 1 ? OCT_DIR_UP : OCT_DIR_FRONT,
                               {n - 4, n - 3},
                               n - 5},
                              {"LEFT", OCT_DIR_LEFT, {2, n - 3}, -1}};
        for (const Case& c : cases) {
            std::fill(Du.begin(), Du.end(), 0.0);
            if (axis == 1)
                ko_deriv21_y(Du.data(), u.data(), h, sz, 1u << c.dir);
            else
                ko_deriv21_z(Du.data(), u.data(), h, sz, 1u << c.dir);
            int miss = 0, checked = 0;
            for (int a = 2; a < n - 2; a++)      // x index
                for (int b = 2; b < n - 2; b++)  // the other transverse index
                    for (int row : c.rows) {
                        const int j  = axis == 1 ? row : b,
                                  k  = axis == 1 ? b : row;
                        const int pp = IDX(a, j, k);
                        const int p0 = pp + (c.base - row) * stride;
                        const double want =
                            c.base >= 0
                                ? (u[p0] - 2.0 * u[p0 + stride] +
                                   u[p0 + 2 * stride]) /
                                      (4.0 * h)
                                : -1.0 / 16.0 / h *
                                      (u[pp - 2 * stride] -
                                       4.0 * u[pp - stride] + 6.0 * u[pp] -
                                       4.0 * u[pp + stride] +
                                       u[pp + 2 * stride]);
                        // same formula as derivs.cpp; 1e-12 allows FMA or
                        // vectorized evaluation there (values are O(10))
                        miss += std::fabs(Du[pp] - want) >
                                1e-12 * (1.0 + std::fabs(want));
                        checked++;
                    }
            std::printf(
                "  ko_deriv21_%c, %s face: %d of %d closure-row values "
                "wrong\n",
                axis == 1 ? 'y' : 'z', c.name, miss, checked);
            bad += miss;
        }
    }
    std::printf("testKODeriv21Faces: %s\n", bad ? "FAIL" : "PASS");
    return bad ? 1 : 0;
}

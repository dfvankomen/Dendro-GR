// Oracle: the closure as coded. On a 13x11x9 random block, ko_deriv21_y/z rows
// 2, 3 (low) or n-4, n-3 (high) equal (u0 - 2 u1 + u2) / (4 h) on the nearest
// three points under the block's own face flag (DOWN/UP, BACK/FRONT); with an x
// flag (LEFT) or none, every row keeps the 5-point interior stencil. Catches:
// flag or axis-stride mix-ups, not closure accuracy. Usage: testKODeriv21Faces

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "derivs.h"

int main() {
    // non-cubic, so a stride taken from the wrong axis fails
    const int N[3] = {13, 11, 9}, nx = N[0], ny = N[1];  // IDX() reads nx, ny
    const unsigned int sz[3] = {13, 11, 9};
    const double h           = 0.1;
    std::vector<double> u(N[0] * N[1] * N[2]), Du(u.size());
    std::mt19937 gen(3);
    std::uniform_real_distribution<double> r(-1.0, 1.0);
    for (double& v : u) v = r(gen);

    int bad = 0;
    for (int axis = 1; axis <= 2; axis++) {
        const int n      = N[axis];
        const int other  = axis == 1 ? 2 : 1;  // the transverse non-x axis
        const int stride = axis == 1 ? nx : nx * ny;
        std::vector<int> all;
        for (int row = 2; row < n - 2; row++) all.push_back(row);
        struct Case {
            const char* name;
            unsigned int flags;
            std::vector<int> rows;
            bool closure;
            int base;  // first point of the closure stencil
        };
        const Case cases[] = {{axis == 1 ? "DOWN" : "BACK",
                               1u << (axis == 1 ? OCT_DIR_DOWN : OCT_DIR_BACK),
                               {2, 3},
                               true,
                               2},
                              {axis == 1 ? "UP" : "FRONT",
                               1u << (axis == 1 ? OCT_DIR_UP : OCT_DIR_FRONT),
                               {n - 4, n - 3},
                               true,
                               n - 5},
                              {"LEFT", 1u << OCT_DIR_LEFT, all, false, 0},
                              {"no", 0u, all, false, 0}};
        for (const Case& c : cases) {
            std::fill(Du.begin(), Du.end(), 0.0);
            if (axis == 1)
                ko_deriv21_y(Du.data(), u.data(), h, sz, c.flags);
            else
                ko_deriv21_z(Du.data(), u.data(), h, sz, c.flags);
            int miss = 0, checked = 0;
            for (int a = 2; a < N[0] - 2; a++)          // x index
                for (int b = 2; b < N[other] - 2; b++)  // transverse index
                    for (int row : c.rows) {
                        const int j  = axis == 1 ? row : b,
                                  k  = axis == 1 ? b : row;
                        const int pp = IDX(a, j, k);
                        const int p0 = pp + (c.base - row) * stride;
                        const double want =
                            c.closure
                                ? (u[p0] - 2.0 * u[p0 + stride] +
                                   u[p0 + 2 * stride]) /
                                      (4.0 * h)
                                : -1.0 / 16.0 / h *
                                      (u[pp - 2 * stride] -
                                       4.0 * u[pp - stride] + 6.0 * u[pp] -
                                       4.0 * u[pp + stride] +
                                       u[pp + 2 * stride]);
                        // same formula as derivs.cpp; 1e-12 allows FMA or
                        // vectorized evaluation there
                        miss += std::fabs(Du[pp] - want) >
                                1e-12 * (1.0 + std::fabs(want));
                        checked++;
                    }
            std::printf("  ko_deriv21_%c, %s face: %d of %d rows wrong\n",
                        axis == 1 ? 'y' : 'z', c.name, miss, checked);
            bad += miss;
        }
    }
    std::printf("testKODeriv21Faces: %s\n", bad ? "FAIL" : "PASS");
    return bad ? 1 : 0;
}

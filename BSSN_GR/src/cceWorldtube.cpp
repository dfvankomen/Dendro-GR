#include "cceWorldtube.h"

#include <algorithm>
#include <cmath>

#include <gsl/gsl_integration.h>

// bssn2adm.h / bssn2adm_derivs.h are pulled from dendrolib/GR/include, which
// this file's own #include below (bssn2adm.h) re-includes as an inline
// fragment inside bssn_point_to_adm_row() -- resolved via
// BSSN_INCLUDE_DIRS in BSSN_GR/CMakeLists.txt (${DENDRO_dendrolib_DIR}/GR/include).

namespace bssn {
namespace cce {

void generate_swsh_angles(const unsigned int l_max, std::vector<double>& theta,
                          std::vector<double>& phi) {
    const unsigned int n_theta = num_theta_points(l_max);
    const unsigned int n_phi   = num_phi_points(l_max);

    // Gauss-Legendre abscissas on [-1, 1] via GSL. We do not rely on any
    // particular ordering GSL returns them in -- we collect all of them,
    // then sort by DESCENDING x (= ascending theta = acos(x)) ourselves,
    // so that theta[0] is nearest the north pole (theta=0) and theta[n-1]
    // is nearest the south pole (theta=pi), matching the ring-index
    // convention documented in cceWorldtube.h / Dendro_CCE_v2.0.md Sec 6.
    gsl_integration_glfixed_table* table =
        gsl_integration_glfixed_table_alloc(n_theta);

    std::vector<double> cos_theta(n_theta);
    for (unsigned int i = 0; i < n_theta; i++) {
        double xi, wi;
        gsl_integration_glfixed_point(-1.0, 1.0, i, &xi, &wi, table);
        cos_theta[i] = xi;
    }
    gsl_integration_glfixed_table_free(table);

    std::sort(cos_theta.begin(), cos_theta.end(), std::greater<double>());

    theta.resize(n_theta);
    for (unsigned int i = 0; i < n_theta; i++) {
        theta[i] = std::acos(cos_theta[i]);
    }

    phi.resize(n_phi);
    for (unsigned int j = 0; j < n_phi; j++) {
        phi[j] = 2.0 * M_PI * static_cast<double>(j) /
                 static_cast<double>(n_phi);
    }
}

void build_worldtube_domain_coords(const std::vector<double>& theta,
                                   const std::vector<double>& phi,
                                   const double radius, const Point& center,
                                   std::vector<double>& domain_coords) {
    const unsigned int n_theta = static_cast<unsigned int>(theta.size());
    const unsigned int n_phi   = static_cast<unsigned int>(phi.size());
    const unsigned int l_max   = n_theta - 1;  // n_theta = l_max + 1

    domain_coords.assign(3u * n_theta * n_phi, 0.0);

    for (unsigned int t = 0; t < n_theta; t++) {
        const double sin_th = std::sin(theta[t]);
        const double cos_th = std::cos(theta[t]);
        for (unsigned int p = 0; p < n_phi; p++) {
            const unsigned int idx = collocation_offset(t, p, l_max);
            const double x =
                center.x() + radius * sin_th * std::cos(phi[p]);
            const double y =
                center.y() + radius * sin_th * std::sin(phi[p]);
            const double z = center.z() + radius * cos_th;
            domain_coords[3 * idx + 0] = x;
            domain_coords[3 * idx + 1] = y;
            domain_coords[3 * idx + 2] = z;
        }
    }
}

const std::vector<std::string>& cce_dataset_names() {
    // clang-format off
    static const std::vector<std::string> names = {
        // g_ij (6)
        "gxx", "gxy", "gxz", "gyy", "gyz", "gzz",
        // d/dx g_ij (6)
        "Dxgxx", "Dxgxy", "Dxgxz", "Dxgyy", "Dxgyz", "Dxgzz",
        // d/dy g_ij (6)
        "Dygxx", "Dygxy", "Dygxz", "Dygyy", "Dygyz", "Dygzz",
        // d/dz g_ij (6)
        "Dzgxx", "Dzgxy", "Dzgxz", "Dzgyy", "Dzgyz", "Dzgzz",
        // Shift^i (3)
        "Shiftx", "Shifty", "Shiftz",
        // d/dx Shift^i (3)
        "DxShiftx", "DxShifty", "DxShiftz",
        // d/dy Shift^i (3)
        "DyShiftx", "DyShifty", "DyShiftz",
        // d/dz Shift^i (3)
        "DzShiftx", "DzShifty", "DzShiftz",
        // Lapse (1)
        "Lapse",
        // d/dx,dy,dz Lapse (3)
        "DxLapse", "DyLapse", "DzLapse",
        // K_ij (6)
        "Kxx", "Kxy", "Kxz", "Kyy", "Kyz", "Kzz",
        // AuxiliaryShift B^i (3)
        "AuxiliaryShiftx", "AuxiliaryShifty", "AuxiliaryShiftz"};
    // clang-format on
    return names;
}

namespace {
// Column index into cce_dataset_names() / row_values for each quantity.
// Must stay in exact sync with cce_dataset_names()'s literal order above.
enum DatasetIndex {
    GXX = 0, GXY, GXZ, GYY, GYZ, GZZ,
    DXGXX, DXGXY, DXGXZ, DXGYY, DXGYZ, DXGZZ,
    DYGXX, DYGXY, DYGXZ, DYGYY, DYGYZ, DYGZZ,
    DZGXX, DZGXY, DZGXZ, DZGYY, DZGYZ, DZGZZ,
    SHIFTX, SHIFTY, SHIFTZ,
    DXSHIFTX, DXSHIFTY, DXSHIFTZ,
    DYSHIFTX, DYSHIFTY, DYSHIFTZ,
    DZSHIFTX, DZSHIFTY, DZSHIFTZ,
    LAPSE,
    DXLAPSE, DYLAPSE, DZLAPSE,
    KXX, KXY, KXZ, KYY, KYZ, KZZ,
    AUXBX, AUXBY, AUXBZ,
    NUM_DATASETS
};

// Unpacks the 6 independent symmetric-tensor components (xx,xy,xz,yy,yz,zz
// order, matching the U_SYMGT0-5 / U_SYMAT0-5 convention confirmed from
// TwoPunctures.cpp's use of adm2bssn.h) into a full 3x3 symmetric matrix.
void unpack_sym(const double c0, const double c1, const double c2,
                const double c3, const double c4, const double c5,
                double m[3][3]) {
    m[0][0] = c0;
    m[0][1] = m[1][0] = c1;
    m[0][2] = m[2][0] = c2;
    m[1][1] = c3;
    m[1][2] = m[2][1] = c4;
    m[2][2] = c5;
}
}  // namespace

void bssn_point_to_adm_row(const double* base_values,
                           const double* deriv_values,
                           const unsigned int point_column,
                           std::vector<std::vector<double>>& row_values) {
    // --- unpack base (0th-derivative) BSSN quantities at this point ---
    const double chi = base_values[CHI];
    const double trK = base_values[TRK];

    double gtd[3][3];
    unpack_sym(base_values[GT0], base_values[GT1], base_values[GT2],
              base_values[GT3], base_values[GT4], base_values[GT5], gtd);

    double Atd[3][3];
    unpack_sym(base_values[AT0], base_values[AT1], base_values[AT2],
              base_values[AT3], base_values[AT4], base_values[AT5], Atd);

    // --- pointwise BSSN -> ADM conversion (produces gd[3][3], Kd[3][3]) ---
#include "bssn2adm.h"

    // --- unpack Cartesian derivatives needed for d_k g_ij ---
    double dchi[3];
    for (int k = 0; k < 3; k++) {
        dchi[k] = deriv_values[deriv_dof(D_CHI, k)];
    }
    double dgtd[3][3][3];
    for (int k = 0; k < 3; k++) {
        const DerivScalar gt_scalars[6] = {D_GT0, D_GT1, D_GT2,
                                           D_GT3, D_GT4, D_GT5};
        double comps[6];
        for (int c = 0; c < 6; c++) {
            comps[c] = deriv_values[deriv_dof(gt_scalars[c],
                                              static_cast<unsigned int>(k))];
        }
        unpack_sym(comps[0], comps[1], comps[2], comps[3], comps[4],
                  comps[5], dgtd[k]);
    }

    // --- chain-rule derivative of the physical metric (produces
    //     dgd[3][3][3]) ---
#include "bssn2adm_derivs.h"

    // --- lapse / shift / auxiliary shift: identical in BSSN and ADM, no
    //     conversion needed, just pass through ---
    const double alpha = base_values[ALPHA];
    const double beta[3] = {base_values[BETA0], base_values[BETA1],
                            base_values[BETA2]};
    const double auxB[3] = {base_values[AUXB0], base_values[AUXB1],
                            base_values[AUXB2]};
    double dalpha[3];
    for (int k = 0; k < 3; k++) {
        dalpha[k] = deriv_values[deriv_dof(D_ALPHA, k)];
    }
    double dbeta[3][3];  // dbeta[k][i] = d_k beta^i
    {
        const DerivScalar beta_scalars[3] = {D_BETA0, D_BETA1, D_BETA2};
        for (int k = 0; k < 3; k++) {
            for (int i = 0; i < 3; i++) {
                dbeta[k][i] = deriv_values[deriv_dof(
                    beta_scalars[i], static_cast<unsigned int>(k))];
            }
        }
    }

    // --- assemble the 49 dataset values at this collocation point ---
    auto set = [&](const DatasetIndex idx, const double value) {
        row_values[idx][point_column] = value;
    };

    set(GXX, gd[0][0]); set(GXY, gd[0][1]); set(GXZ, gd[0][2]);
    set(GYY, gd[1][1]); set(GYZ, gd[1][2]); set(GZZ, gd[2][2]);

    set(DXGXX, dgd[0][0][0]); set(DXGXY, dgd[0][0][1]);
    set(DXGXZ, dgd[0][0][2]); set(DXGYY, dgd[0][1][1]);
    set(DXGYZ, dgd[0][1][2]); set(DXGZZ, dgd[0][2][2]);

    set(DYGXX, dgd[1][0][0]); set(DYGXY, dgd[1][0][1]);
    set(DYGXZ, dgd[1][0][2]); set(DYGYY, dgd[1][1][1]);
    set(DYGYZ, dgd[1][1][2]); set(DYGZZ, dgd[1][2][2]);

    set(DZGXX, dgd[2][0][0]); set(DZGXY, dgd[2][0][1]);
    set(DZGXZ, dgd[2][0][2]); set(DZGYY, dgd[2][1][1]);
    set(DZGYZ, dgd[2][1][2]); set(DZGZZ, dgd[2][2][2]);

    set(SHIFTX, beta[0]); set(SHIFTY, beta[1]); set(SHIFTZ, beta[2]);

    set(DXSHIFTX, dbeta[0][0]); set(DXSHIFTY, dbeta[0][1]);
    set(DXSHIFTZ, dbeta[0][2]);
    set(DYSHIFTX, dbeta[1][0]); set(DYSHIFTY, dbeta[1][1]);
    set(DYSHIFTZ, dbeta[1][2]);
    set(DZSHIFTX, dbeta[2][0]); set(DZSHIFTY, dbeta[2][1]);
    set(DZSHIFTZ, dbeta[2][2]);

    set(LAPSE, alpha);
    set(DXLAPSE, dalpha[0]); set(DYLAPSE, dalpha[1]); set(DZLAPSE, dalpha[2]);

    set(KXX, Kd[0][0]); set(KXY, Kd[0][1]); set(KXZ, Kd[0][2]);
    set(KYY, Kd[1][1]); set(KYZ, Kd[1][2]); set(KZZ, Kd[2][2]);

    set(AUXBX, auxB[0]); set(AUXBY, auxB[1]); set(AUXBZ, auxB[2]);
}

}  // namespace cce
}  // namespace bssn

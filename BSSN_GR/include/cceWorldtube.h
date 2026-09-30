/**
 * @file cceWorldtube.h
 * @brief Pure-math / pure-IO helpers for writing a Cauchy-Characteristic
 *        Extraction (CCE) worldtube HDF5 file, consumable by SpECTRE's
 *        PreprocessCceWorldtube (AdmMetricNodal format, first-order shift
 *        driver). See Dendro_CCE_v2.0.md for the full design and the exact
 *        data contract this implements, with citations.
 *
 * This header intentionally contains NO ot::Mesh / BSSNCtx dependent code
 * (no unzip/zip, no interpolateToCoords calls) -- that orchestration lives
 * in BSSNCtx::writeCceWorldtube() (bssnCtx.h/.cpp), because it needs
 * BSSNCtx's protected zip()/unzip() methods. Everything here is testable in
 * isolation: SWSH grid generation, the Cartesian query-point construction,
 * and the per-point BSSN->ADM conversion + HDF5 row assembly.
 */

#ifndef BSSN_CCE_WORLDTUBE_H
#define BSSN_CCE_WORLDTUBE_H

#include <string>
#include <vector>

#include "cceH5Dat.h"
#include "grDef.h"
#include "point.h"

namespace bssn {
namespace cce {

// ---------------------------------------------------------------------
// SWSH collocation grid (matches SpECTRE's libsharp Gauss-Legendre grid,
// src/NumericalAlgorithms/SpinWeightedSphericalHarmonics/SwshCollocation.hpp
// in the spectre repo -- see Dendro_CCE_v2.0.md Section 6 for full
// citations and the UNVERIFIED flag on the phi-fastest/theta-slowest
// ordering assumed here).
// ---------------------------------------------------------------------

inline unsigned int num_theta_points(const unsigned int l_max) {
    return l_max + 1;
}
inline unsigned int num_phi_points(const unsigned int l_max) {
    return 2 * l_max + 1;
}
inline unsigned int num_collocation_points(const unsigned int l_max) {
    return num_theta_points(l_max) * num_phi_points(l_max);
}
// offset = phi_index + num_phi_points(l_max) * theta_ring_index
// (phi fastest / inner loop, theta slowest / outer loop).
inline unsigned int collocation_offset(const unsigned int theta_ring,
                                       const unsigned int phi_index,
                                       const unsigned int l_max) {
    return phi_index + num_phi_points(l_max) * theta_ring;
}

/**
 * @brief Fills `theta` (size num_theta_points(l_max)) with the Gauss-
 *        Legendre collocation angles in [0, pi], ring 0 nearest theta=0,
 *        ascending; and `phi` (size num_phi_points(l_max)) with the
 *        equally-spaced angles phi_j = 2*pi*j/(2*l_max+1), j=0,...,2*l_max
 *        (phi_0 = 0, no phase offset). Uses GSL's
 *        gsl_integration_glfixed_table for the Gauss-Legendre abscissas,
 *        then sorts them into descending cos(theta) (= ascending theta)
 *        order itself, rather than relying on any particular ordering
 *        convention GSL returns them in.
 */
void generate_swsh_angles(unsigned int l_max, std::vector<double>& theta,
                          std::vector<double>& phi);

/**
 * @brief Builds the flat Cartesian (x,y,z) query-point list for
 *        ot::da::interpolateToCoords, for a sphere of coordinate `radius`
 *        centered at `center`, on the (theta, phi) grid from
 *        generate_swsh_angles(). `domain_coords` is resized to
 *        3*theta.size()*phi.size(), laid out as
 *        domain_coords[3*collocation_offset(...) + {0,1,2}] = {x,y,z},
 *        mirroring AEH_BHaHAHA::fill_domain_coords's layout convention
 *        (dendrolib/src/aeh_bhahaha.cpp).
 */
void build_worldtube_domain_coords(const std::vector<double>& theta,
                                   const std::vector<double>& phi,
                                   double radius, const Point& center,
                                   std::vector<double>& domain_coords);

// ---------------------------------------------------------------------
// Field bookkeeping: which BSSN evolved variables get interpolated
// directly (0th derivative), and which additionally need volume Cartesian
// derivatives computed before interpolation. Order fixes the dof indexing
// used by BSSNCtx::writeCceWorldtube() when it calls
// ot::da::interpolateToCoords once per field, matching the per-field-loop
// pattern in AEH_BHaHAHA::interpolate_metric_data.
// ---------------------------------------------------------------------

// clang-format off
enum BaseField {
    CHI = 0, GT0, GT1, GT2, GT3, GT4, GT5, AT0, AT1, AT2, AT3, AT4, AT5,
    TRK, ALPHA, BETA0, BETA1, BETA2, AUXB0, AUXB1, AUXB2, NUM_BASE_FIELDS
};
// clang-format on

// VAR indices, in BaseField order -- passed to ot::da::interpolateToCoords
// as varData[CCE_BASE_INDICES[f]].
const std::vector<int> CCE_BASE_INDICES = {
    bssn::VAR::U_CHI,    bssn::VAR::U_SYMGT0, bssn::VAR::U_SYMGT1,
    bssn::VAR::U_SYMGT2, bssn::VAR::U_SYMGT3, bssn::VAR::U_SYMGT4,
    bssn::VAR::U_SYMGT5, bssn::VAR::U_SYMAT0, bssn::VAR::U_SYMAT1,
    bssn::VAR::U_SYMAT2, bssn::VAR::U_SYMAT3, bssn::VAR::U_SYMAT4,
    bssn::VAR::U_SYMAT5, bssn::VAR::U_K,      bssn::VAR::U_ALPHA,
    bssn::VAR::U_BETA0,  bssn::VAR::U_BETA1,  bssn::VAR::U_BETA2,
    bssn::VAR::U_B0,     bssn::VAR::U_B1,     bssn::VAR::U_B2};

// clang-format off
// Scalars needing volume Cartesian derivatives d/dx,d/dy,d/dz BEFORE
// interpolation (chi, the 6 independent gtd components, alpha, the 3 beta
// components -- 11 scalars x 3 directions = 33 derivative dof slots).
// AdmMetricNodal does NOT need derivatives of Atd, trK, or B^i (see
// Dendro_CCE_v2.0.md Section 6.1's schema table), so those are excluded.
enum DerivScalar {
    D_CHI = 0, D_GT0, D_GT1, D_GT2, D_GT3, D_GT4, D_GT5, D_ALPHA, D_BETA0,
    D_BETA1, D_BETA2, NUM_DERIV_SCALARS
};
// clang-format on
const unsigned int NUM_DERIV_DOFS = 3u * NUM_DERIV_SCALARS;  // = 33

// VAR index of the field each DerivScalar differentiates.
const std::vector<int> DERIV_SCALAR_SOURCE_VAR = {
    bssn::VAR::U_CHI,   bssn::VAR::U_SYMGT0, bssn::VAR::U_SYMGT1,
    bssn::VAR::U_SYMGT2, bssn::VAR::U_SYMGT3, bssn::VAR::U_SYMGT4,
    bssn::VAR::U_SYMGT5, bssn::VAR::U_ALPHA, bssn::VAR::U_BETA0,
    bssn::VAR::U_BETA1, bssn::VAR::U_BETA2};

// derivative dof index for (scalar, direction), direction: 0=x,1=y,2=z
inline unsigned int deriv_dof(const DerivScalar scalar,
                              const unsigned int direction) {
    return 3u * static_cast<unsigned int>(scalar) + direction;
}

/**
 * @brief The 49 AdmMetricNodal (first-order driver) HDF5 dataset base
 *        names, in the fixed order BSSNCtx::writeCceWorldtube() assembles
 *        columns in -- see Dendro_CCE_v2.0.md Section 6.1 for the schema
 *        this implements (g/Dxg/Dyg/Dzg: 6 each = 24; Shift/DxShift/
 *        DyShift/DzShift: 3 each = 12; Lapse: 1, DxLapse/DyLapse/DzLapse:
 *        3 = 4; K: 6; AuxiliaryShift: 3. Total 49).
 */
const std::vector<std::string>& cce_dataset_names();

/**
 * @brief Given the interpolated BSSN-frame values at ONE collocation point
 *        (base fields, in BaseField order, and derivative fields, in
 *        deriv_dof() order), computes the physical ADM quantities via
 *        bssn2adm.h / bssn2adm_derivs.h and writes them into `row_values`
 *        (size 49, matching cce_dataset_names()' order) at column index
 *        `point_column` (i.e. row_values[d][point_column] = value of
 *        dataset d at this point -- row_values is indexed [dataset][col],
 *        caller owns allocation so this can be called once per point
 *        without re-allocating each time).
 */
void bssn_point_to_adm_row(const double* base_values,
                           const double* deriv_values,
                           unsigned int point_column,
                           std::vector<std::vector<double>>& row_values);

}  // namespace cce
}  // namespace bssn

#endif  // BSSN_CCE_WORLDTUBE_H

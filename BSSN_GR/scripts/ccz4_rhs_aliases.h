const double p_expo = 1.0;
const bool ccz4_diag_print = false;
auto ccz4_diag_rank = []() -> int { return 0; };

const double *const psi = chi;
double *const psi_rhs   = chi_rhs;

const double *const gammat00 = gt0;
const double *const gammat01 = gt1;
const double *const gammat02 = gt2;
const double *const gammat11 = gt3;
const double *const gammat12 = gt4;
const double *const gammat22 = gt5;

const double *const At00 = At0;
const double *const At01 = At1;
const double *const At02 = At2;
const double *const At11 = At3;
const double *const At12 = At4;
const double *const At22 = At5;

double *const alpha_rhs = a_rhs;
double *const beta0_rhs = b_rhs0;
double *const beta1_rhs = b_rhs1;
double *const beta2_rhs = b_rhs2;

double *const gammat00_rhs = gt_rhs00;
double *const gammat01_rhs = gt_rhs01;
double *const gammat02_rhs = gt_rhs02;
double *const gammat11_rhs = gt_rhs11;
double *const gammat12_rhs = gt_rhs12;
double *const gammat22_rhs = gt_rhs22;

double *const At00_rhs = At_rhs00;
double *const At01_rhs = At_rhs01;
double *const At02_rhs = At_rhs02;
double *const At11_rhs = At_rhs11;
double *const At12_rhs = At_rhs12;
double *const At22_rhs = At_rhs22;

double *const B0_rhs = B_rhs0;
double *const B1_rhs = B_rhs1;
double *const B2_rhs = B_rhs2;

#define CCZ4_ALIAS_DERIVS(name, base)                                           \
    double *const grad_0_##name    = grad_0_##base;                             \
    double *const grad_1_##name    = grad_1_##base;                             \
    double *const grad_2_##name    = grad_2_##base;                             \
    double *const grad2_0_0_##name = grad2_0_0_##base;                          \
    double *const grad2_0_1_##name = grad2_0_1_##base;                          \
    double *const grad2_0_2_##name = grad2_0_2_##base;                          \
    double *const grad2_1_1_##name = grad2_1_1_##base;                          \
    double *const grad2_1_2_##name = grad2_1_2_##base;                          \
    double *const grad2_2_2_##name = grad2_2_2_##base

CCZ4_ALIAS_DERIVS(psi, chi);
CCZ4_ALIAS_DERIVS(gammat00, gt0);
CCZ4_ALIAS_DERIVS(gammat01, gt1);
CCZ4_ALIAS_DERIVS(gammat02, gt2);
CCZ4_ALIAS_DERIVS(gammat11, gt3);
CCZ4_ALIAS_DERIVS(gammat12, gt4);
CCZ4_ALIAS_DERIVS(gammat22, gt5);
#undef CCZ4_ALIAS_DERIVS

#define CCZ4_ALIAS_FIRST_DERIVS(name, base)                                     \
    double *const grad_0_##name = grad_0_##base;                                \
    double *const grad_1_##name = grad_1_##base;                                \
    double *const grad_2_##name = grad_2_##base

CCZ4_ALIAS_FIRST_DERIVS(At00, At0);
CCZ4_ALIAS_FIRST_DERIVS(At01, At1);
CCZ4_ALIAS_FIRST_DERIVS(At02, At2);
CCZ4_ALIAS_FIRST_DERIVS(At11, At3);
CCZ4_ALIAS_FIRST_DERIVS(At12, At4);
CCZ4_ALIAS_FIRST_DERIVS(At22, At5);
#undef CCZ4_ALIAS_FIRST_DERIVS

#ifndef BSSN_USE_ADVECTIVE_DERIVS
double *const agrad_0_beta0    = grad_0_beta0;
double *const agrad_1_beta0    = grad_1_beta0;
double *const agrad_2_beta0    = grad_2_beta0;
double *const agrad_0_beta1    = grad_0_beta1;
double *const agrad_1_beta1    = grad_1_beta1;
double *const agrad_2_beta1    = grad_2_beta1;
double *const agrad_0_beta2    = grad_0_beta2;
double *const agrad_1_beta2    = grad_1_beta2;
double *const agrad_2_beta2    = grad_2_beta2;
double *const agrad_0_B0        = grad_0_B0;
double *const agrad_1_B0        = grad_1_B0;
double *const agrad_2_B0        = grad_2_B0;
double *const agrad_0_B1        = grad_0_B1;
double *const agrad_1_B1        = grad_1_B1;
double *const agrad_2_B1        = grad_2_B1;
double *const agrad_0_B2        = grad_0_B2;
double *const agrad_1_B2        = grad_1_B2;
double *const agrad_2_B2        = grad_2_B2;
double *const agrad_0_Gammahat0 = grad_0_Gammahat0;
double *const agrad_1_Gammahat0 = grad_1_Gammahat0;
double *const agrad_2_Gammahat0 = grad_2_Gammahat0;
double *const agrad_0_Gammahat1 = grad_0_Gammahat1;
double *const agrad_1_Gammahat1 = grad_1_Gammahat1;
double *const agrad_2_Gammahat1 = grad_2_Gammahat1;
double *const agrad_0_Gammahat2 = grad_0_Gammahat2;
double *const agrad_1_Gammahat2 = grad_1_Gammahat2;
double *const agrad_2_Gammahat2 = grad_2_Gammahat2;
#endif

double *const ccz4_workspace = __mem_pool->allocate(50 * n);
double *const Z0                         = ccz4_workspace + 0 * n;
double *const Z1                         = ccz4_workspace + 1 * n;
double *const Z2                         = ccz4_workspace + 2 * n;
double *const ccz4_D0310                 = ccz4_workspace + 3 * n;
double *const ccz4_D0312                 = ccz4_workspace + 4 * n;
double *const ccz4_D0314                 = ccz4_workspace + 5 * n;
double *const ccz4_D0293                 = ccz4_workspace + 6 * n;
double *const ccz4_D0426                 = ccz4_workspace + 7 * n;
double *const ccz4_D0521                 = ccz4_workspace + 8 * n;
double *const ccz4_D0585                 = ccz4_workspace + 9 * n;
double *const ccz4_D0931                 = ccz4_workspace + 10 * n;
double *const ccz4_D0521_norm            = ccz4_workspace + 11 * n;
double *const ccz4_mD0293_norm           = ccz4_workspace + 12 * n;
double *const ccz4_mD0930_norm           = ccz4_workspace + 13 * n;
double *const ccz4_D0585_norm            = ccz4_workspace + 14 * n;
double *const ccz4_D0426_norm            = ccz4_workspace + 15 * n;
double *const grad_0_ccz4_D0310          = ccz4_workspace + 16 * n;
double *const grad_1_ccz4_D0310          = ccz4_workspace + 17 * n;
double *const grad_2_ccz4_D0310          = ccz4_workspace + 18 * n;
double *const grad_0_ccz4_D0312          = ccz4_workspace + 19 * n;
double *const grad_1_ccz4_D0312          = ccz4_workspace + 20 * n;
double *const grad_2_ccz4_D0312          = ccz4_workspace + 21 * n;
double *const grad_0_ccz4_D0314          = ccz4_workspace + 22 * n;
double *const grad_1_ccz4_D0314          = ccz4_workspace + 23 * n;
double *const grad_2_ccz4_D0314          = ccz4_workspace + 24 * n;
double *const grad_0_ccz4_D0293          = ccz4_workspace + 25 * n;
double *const grad_0_ccz4_D0426          = ccz4_workspace + 26 * n;
double *const grad_1_ccz4_D0426          = ccz4_workspace + 27 * n;
double *const grad_2_ccz4_D0426          = ccz4_workspace + 28 * n;
double *const grad_0_ccz4_D0521          = ccz4_workspace + 29 * n;
double *const grad_1_ccz4_D0521          = ccz4_workspace + 30 * n;
double *const grad_2_ccz4_D0521          = ccz4_workspace + 31 * n;
double *const grad_0_ccz4_D0585          = ccz4_workspace + 32 * n;
double *const grad_1_ccz4_D0585          = ccz4_workspace + 33 * n;
double *const grad_2_ccz4_D0585          = ccz4_workspace + 34 * n;
double *const grad_0_ccz4_D0931          = ccz4_workspace + 35 * n;
double *const grad_1_ccz4_D0931          = ccz4_workspace + 36 * n;
double *const grad_2_ccz4_D0931          = ccz4_workspace + 37 * n;
double *const grad_1_ccz4_D0521_norm     = ccz4_workspace + 38 * n;
double *const grad_0_ccz4_mD0293_norm    = ccz4_workspace + 39 * n;
double *const grad_2_ccz4_mD0930_norm    = ccz4_workspace + 40 * n;
double *const grad_0_ccz4_D0585_norm     = ccz4_workspace + 41 * n;
double *const grad_2_ccz4_D0426_norm     = ccz4_workspace + 42 * n;

unsigned int ccz4_invalid_padded_points = 0;

for (unsigned int qq = 0; qq < n; qq++) {
    const double g[3][3] = {{gammat00[qq], gammat01[qq], gammat02[qq]},
                            {gammat01[qq], gammat11[qq], gammat12[qq]},
                            {gammat02[qq], gammat12[qq], gammat22[qq]}};
    const double det = g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[1][2]) -
                       g[0][1] * (g[0][1] * g[2][2] - g[0][2] * g[1][2]) +
                       g[0][2] * (g[0][1] * g[1][2] - g[0][2] * g[1][1]);

    const unsigned int ccz4_i = qq % nx;
    const unsigned int ccz4_j = (qq / nx) % ny;
    const unsigned int ccz4_k = qq / (nx * ny);
    const bool ccz4_is_rhs_point =
        (ccz4_i >= PW && ccz4_i < nx - PW && ccz4_j >= PW &&
         ccz4_j < ny - PW && ccz4_k >= PW && ccz4_k < nz - PW);
    const bool ccz4_invalid_metric_or_psi =
        (!std::isfinite(det) || fabs(det) < 1.0e-300 ||
         !std::isfinite(psi[qq]) || psi[qq] <= 0.0);

    // Dendro unzip padding can contain zero or otherwise invalid temporary
    // data outside the physical RHS loop. These CCZ4 wrapper temps require an
    // inverse metric and powers of psi, so padded invalid points must be made
    // finite before derivative/RHS diagnostics scan full block arrays. This
    // branch is fatal on valid physical RHS points and does not change the
    // generated CCZ4 equations on points used by the RHS loop.
    if (ccz4_invalid_metric_or_psi) {
        if (ccz4_is_rhs_point) {
            std::cerr << "[CCZ4 RHS DIAG] invalid physical CCZ4 temp input"
                      << " pp=" << qq << " ijk=(" << ccz4_i << ","
                      << ccz4_j << "," << ccz4_k << ") det=" << det
                      << " psi=" << psi[qq] << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 913);
        }
        ccz4_invalid_padded_points++;
        const double safe_Gammahat0 =
            std::isfinite(Gammahat0[qq]) ? Gammahat0[qq] : 0.0;
        const double safe_Gammahat1 =
            std::isfinite(Gammahat1[qq]) ? Gammahat1[qq] : 0.0;
        const double safe_Gammahat2 =
            std::isfinite(Gammahat2[qq]) ? Gammahat2[qq] : 0.0;
        Z0[qq] = 0.0;
        Z1[qq] = 0.0;
        Z2[qq] = 0.0;
        ccz4_D0310[qq]       = safe_Gammahat0;
        ccz4_D0312[qq]       = safe_Gammahat1;
        ccz4_D0314[qq]       = safe_Gammahat2;
        ccz4_D0293[qq]       = 0.0;
        ccz4_D0426[qq]       = 0.0;
        ccz4_D0521[qq]       = 0.0;
        ccz4_D0585[qq]       = 0.0;
        ccz4_D0931[qq]       = 0.0;
        ccz4_D0521_norm[qq]  = 0.0;
        ccz4_mD0293_norm[qq] = 0.0;
        ccz4_mD0930_norm[qq] = 0.0;
        ccz4_D0585_norm[qq]  = 0.0;
        ccz4_D0426_norm[qq]  = 0.0;
        continue;
    }

    const double invg[3][3] = {
        {(g[1][1] * g[2][2] - g[1][2] * g[1][2]) / det,
         (g[0][2] * g[1][2] - g[0][1] * g[2][2]) / det,
         (g[0][1] * g[1][2] - g[0][2] * g[1][1]) / det},
        {(g[0][2] * g[1][2] - g[0][1] * g[2][2]) / det,
         (g[0][0] * g[2][2] - g[0][2] * g[0][2]) / det,
         (g[0][1] * g[0][2] - g[0][0] * g[1][2]) / det},
        {(g[0][1] * g[1][2] - g[0][2] * g[1][1]) / det,
         (g[0][1] * g[0][2] - g[0][0] * g[1][2]) / det,
         (g[0][0] * g[1][1] - g[0][1] * g[0][1]) / det}};

    const double dg[3][3][3] = {
        {{grad_0_gammat00[qq], grad_0_gammat01[qq], grad_0_gammat02[qq]},
         {grad_0_gammat01[qq], grad_0_gammat11[qq], grad_0_gammat12[qq]},
         {grad_0_gammat02[qq], grad_0_gammat12[qq], grad_0_gammat22[qq]}},
        {{grad_1_gammat00[qq], grad_1_gammat01[qq], grad_1_gammat02[qq]},
         {grad_1_gammat01[qq], grad_1_gammat11[qq], grad_1_gammat12[qq]},
         {grad_1_gammat02[qq], grad_1_gammat12[qq], grad_1_gammat22[qq]}},
        {{grad_2_gammat00[qq], grad_2_gammat01[qq], grad_2_gammat02[qq]},
         {grad_2_gammat01[qq], grad_2_gammat11[qq], grad_2_gammat12[qq]},
         {grad_2_gammat02[qq], grad_2_gammat12[qq], grad_2_gammat22[qq]}}};

    double contracted[3] = {0.0, 0.0, 0.0};
    for (unsigned int u = 0; u < 3; u++) {
        for (unsigned int j = 0; j < 3; j++) {
            for (unsigned int k = 0; k < 3; k++) {
                for (unsigned int l = 0; l < 3; l++) {
                    contracted[u] +=
                        0.5 * invg[j][k] * invg[u][l] *
                        (dg[j][l][k] + dg[k][l][j] - dg[l][j][k]);
                }
            }
        }
    }

    const double diff[3] = {Gammahat0[qq] - contracted[0],
                            Gammahat1[qq] - contracted[1],
                            Gammahat2[qq] - contracted[2]};
    const double z_scale = 0.5 / pow(psi[qq], p_expo);
    Z0[qq] = z_scale * (g[0][0] * diff[0] + g[0][1] * diff[1] +
                         g[0][2] * diff[2]);
    Z1[qq] = z_scale * (g[1][0] * diff[0] + g[1][1] * diff[1] +
                         g[1][2] * diff[2]);
    Z2[qq] = z_scale * (g[2][0] * diff[0] + g[2][1] * diff[1] +
                         g[2][2] * diff[2]);

    const double D0011 = gammat01[qq] * gammat01[qq];
    const double D0012 = gammat00[qq] * gammat11[qq];
    const double D0013 = -D0011 + D0012;
    const double D0014 = gammat12[qq] * gammat12[qq];
    const double D0015 = gammat02[qq] * gammat02[qq];
    const double D0016 = gammat01[qq] * gammat02[qq];
    const double D0017 = D0011 * gammat22[qq] - D0012 * gammat22[qq] +
                         D0014 * gammat00[qq] + D0015 * gammat11[qq] -
                         2.0 * D0016 * gammat12[qq];
    const double D0018 = 1.0 / D0017;
    const double D0020 = -D0015 + gammat00[qq] * gammat22[qq];
    const double D0022 = gammat01[qq] * gammat12[qq] -
                         gammat02[qq] * gammat11[qq];
    const double D0024 = -D0016 + gammat00[qq] * gammat12[qq];
    const double D0025 = gammat01[qq] * gammat22[qq] -
                         gammat02[qq] * gammat12[qq];
    const double D0026 = -D0014 + gammat11[qq] * gammat22[qq];
    const double D0061 = D0022 * Z2[qq] - D0025 * Z1[qq] + D0026 * Z0[qq];
    const double D0062 = -D0061;
    const double D0194 = -D0020 * Z1[qq] + D0024 * Z2[qq] + D0025 * Z0[qq];
    const double D0225 = D0013 * Z2[qq] + D0022 * Z0[qq] - D0024 * Z1[qq];
    const double D0226 = -D0225;
    const double D0291 = pow(psi[qq], 2.0 * p_expo);
    const double D0292 = D0018 * D0291;
    const double D0309 = D0018 * 2.0 * pow(psi[qq], p_expo);
    const double D0930 = D0225 * D0292;

    ccz4_D0310[qq] = D0061 * D0309 + Gammahat0[qq];
    ccz4_D0312[qq] = -D0194 * D0309 + Gammahat1[qq];
    ccz4_D0314[qq] = D0225 * D0309 + Gammahat2[qq];
    ccz4_D0293[qq] = D0061 * D0292;
    ccz4_D0426[qq] = D0226 * D0292;
    ccz4_D0521[qq] = D0194 * D0292;
    ccz4_D0585[qq] = D0062 * D0292;
    ccz4_D0931[qq] = -D0930;

    const double D1021 =
        sqrt(fabs(gammat00[qq]) * fabs(gammat00[qq]) +
             2.0 * fabs(gammat01[qq]) * fabs(gammat01[qq]) +
             2.0 * fabs(gammat02[qq]) * fabs(gammat02[qq]) +
             fabs(gammat11[qq]) * fabs(gammat11[qq]) +
             2.0 * fabs(gammat12[qq]) * fabs(gammat12[qq]) +
             fabs(gammat22[qq]) * fabs(gammat22[qq]));
    ccz4_D0521_norm[qq]  = ccz4_D0521[qq] * D1021;
    ccz4_mD0293_norm[qq] = -ccz4_D0293[qq] * D1021;
    ccz4_mD0930_norm[qq] = -D0930 * D1021;
    ccz4_D0585_norm[qq]  = ccz4_D0585[qq] * D1021;
    ccz4_D0426_norm[qq]  = ccz4_D0426[qq] * D1021;

    const bool ccz4_temp_finite =
        std::isfinite(Z0[qq]) && std::isfinite(Z1[qq]) &&
        std::isfinite(Z2[qq]) && std::isfinite(ccz4_D0310[qq]) &&
        std::isfinite(ccz4_D0312[qq]) && std::isfinite(ccz4_D0314[qq]) &&
        std::isfinite(ccz4_D0293[qq]) && std::isfinite(ccz4_D0426[qq]) &&
        std::isfinite(ccz4_D0521[qq]) && std::isfinite(ccz4_D0585[qq]) &&
        std::isfinite(ccz4_D0931[qq]) &&
        std::isfinite(ccz4_D0521_norm[qq]) &&
        std::isfinite(ccz4_mD0293_norm[qq]) &&
        std::isfinite(ccz4_mD0930_norm[qq]) &&
        std::isfinite(ccz4_D0585_norm[qq]) &&
        std::isfinite(ccz4_D0426_norm[qq]);
    if (!ccz4_temp_finite) {
        if (ccz4_is_rhs_point) {
            std::cerr << "[CCZ4 RHS DIAG] non-finite physical CCZ4 wrapper temp"
                      << " pp=" << qq << " ijk=(" << ccz4_i << ","
                      << ccz4_j << "," << ccz4_k << ") det=" << det
                      << " psi=" << psi[qq] << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 914);
        }
        ccz4_invalid_padded_points++;
        Z0[qq]                 = 0.0;
        Z1[qq]                 = 0.0;
        Z2[qq]                 = 0.0;
        ccz4_D0310[qq]         = 0.0;
        ccz4_D0312[qq]         = 0.0;
        ccz4_D0314[qq]         = 0.0;
        ccz4_D0293[qq]         = 0.0;
        ccz4_D0426[qq]         = 0.0;
        ccz4_D0521[qq]         = 0.0;
        ccz4_D0585[qq]         = 0.0;
        ccz4_D0931[qq]         = 0.0;
        ccz4_D0521_norm[qq]    = 0.0;
        ccz4_mD0293_norm[qq]   = 0.0;
        ccz4_mD0930_norm[qq]   = 0.0;
        ccz4_D0585_norm[qq]    = 0.0;
        ccz4_D0426_norm[qq]    = 0.0;
    }
}

if (ccz4_diag_print && ccz4_invalid_padded_points > 0 &&
    ccz4_diag_rank() == 0) {
    std::cout << "[CCZ4 RHS DIAG] skipped invalid padded/temp points="
              << ccz4_invalid_padded_points
              << " before CCZ4 temp derivative construction" << std::endl;
}

deriv_x(grad_0_ccz4_D0310, ccz4_D0310, hx, sz, bflag);
deriv_y(grad_1_ccz4_D0310, ccz4_D0310, hy, sz, bflag);
deriv_z(grad_2_ccz4_D0310, ccz4_D0310, hz, sz, bflag);
deriv_x(grad_0_ccz4_D0312, ccz4_D0312, hx, sz, bflag);
deriv_y(grad_1_ccz4_D0312, ccz4_D0312, hy, sz, bflag);
deriv_z(grad_2_ccz4_D0312, ccz4_D0312, hz, sz, bflag);
deriv_x(grad_0_ccz4_D0314, ccz4_D0314, hx, sz, bflag);
deriv_y(grad_1_ccz4_D0314, ccz4_D0314, hy, sz, bflag);
deriv_z(grad_2_ccz4_D0314, ccz4_D0314, hz, sz, bflag);
deriv_x(grad_0_ccz4_D0293, ccz4_D0293, hx, sz, bflag);
deriv_x(grad_0_ccz4_D0426, ccz4_D0426, hx, sz, bflag);
deriv_y(grad_1_ccz4_D0426, ccz4_D0426, hy, sz, bflag);
deriv_z(grad_2_ccz4_D0426, ccz4_D0426, hz, sz, bflag);
deriv_x(grad_0_ccz4_D0521, ccz4_D0521, hx, sz, bflag);
deriv_y(grad_1_ccz4_D0521, ccz4_D0521, hy, sz, bflag);
deriv_z(grad_2_ccz4_D0521, ccz4_D0521, hz, sz, bflag);
deriv_x(grad_0_ccz4_D0585, ccz4_D0585, hx, sz, bflag);
deriv_y(grad_1_ccz4_D0585, ccz4_D0585, hy, sz, bflag);
deriv_z(grad_2_ccz4_D0585, ccz4_D0585, hz, sz, bflag);
deriv_x(grad_0_ccz4_D0931, ccz4_D0931, hx, sz, bflag);
deriv_y(grad_1_ccz4_D0931, ccz4_D0931, hy, sz, bflag);
deriv_z(grad_2_ccz4_D0931, ccz4_D0931, hz, sz, bflag);
deriv_y(grad_1_ccz4_D0521_norm, ccz4_D0521_norm, hy, sz, bflag);
deriv_x(grad_0_ccz4_mD0293_norm, ccz4_mD0293_norm, hx, sz, bflag);
deriv_z(grad_2_ccz4_mD0930_norm, ccz4_mD0930_norm, hz, sz, bflag);
deriv_x(grad_0_ccz4_D0585_norm, ccz4_D0585_norm, hx, sz, bflag);
deriv_z(grad_2_ccz4_D0426_norm, ccz4_D0426_norm, hz, sz, bflag);

auto ccz4_grad = [&](const int dir, const char *expr,
                     const unsigned int pp) -> double {
    if (std::strcmp(expr, "DENDRO_0056") == 0) {
        if (dir == 0) return grad_0_alpha[pp] / alpha[pp];
        if (dir == 1) return grad_1_alpha[pp] / alpha[pp];
        if (dir == 2) return grad_2_alpha[pp] / alpha[pp];
    }
    if (std::strcmp(expr, "DENDRO_0064") == 0) {
        if (dir == 0) return grad_0_K[pp] + 2.0 * kappa[4] * grad_0_Theta[pp];
        if (dir == 1) return grad_1_K[pp] + 2.0 * kappa[4] * grad_1_Theta[pp];
        if (dir == 2) return grad_2_K[pp] + 2.0 * kappa[4] * grad_2_Theta[pp];
    }
    if (std::strcmp(expr, "DENDRO_0133") == 0) {
        if (dir == 0) return grad_0_psi[pp] / psi[pp];
        if (dir == 1) return grad_1_psi[pp] / psi[pp];
        if (dir == 2) return grad_2_psi[pp] / psi[pp];
    }
    if (std::strcmp(expr, "DENDRO_0310") == 0 ||
        std::strcmp(expr, "DENDRO_0435") == 0) {
        if (dir == 0) return grad_0_ccz4_D0310[pp];
        if (dir == 1) return grad_1_ccz4_D0310[pp];
        if (dir == 2) return grad_2_ccz4_D0310[pp];
    }
    if (std::strcmp(expr, "DENDRO_0312") == 0) {
        if (dir == 0) return grad_0_ccz4_D0312[pp];
        if (dir == 1) return grad_1_ccz4_D0312[pp];
        if (dir == 2) return grad_2_ccz4_D0312[pp];
    }
    if (std::strcmp(expr, "DENDRO_0314") == 0 ||
        std::strcmp(expr, "DENDRO_0438") == 0) {
        if (dir == 0) return grad_0_ccz4_D0314[pp];
        if (dir == 1) return grad_1_ccz4_D0314[pp];
        if (dir == 2) return grad_2_ccz4_D0314[pp];
    }
    if (std::strcmp(expr, "DENDRO_0426") == 0) {
        if (dir == 0) return grad_0_ccz4_D0426[pp];
        if (dir == 1) return grad_1_ccz4_D0426[pp];
        if (dir == 2) return grad_2_ccz4_D0426[pp];
    }
    if (std::strcmp(expr, "DENDRO_0521") == 0) {
        if (dir == 0) return grad_0_ccz4_D0521[pp];
        if (dir == 1) return grad_1_ccz4_D0521[pp];
        if (dir == 2) return grad_2_ccz4_D0521[pp];
    }
    if (std::strcmp(expr, "DENDRO_0585") == 0) {
        if (dir == 0) return grad_0_ccz4_D0585[pp];
        if (dir == 1) return grad_1_ccz4_D0585[pp];
        if (dir == 2) return grad_2_ccz4_D0585[pp];
    }
    if (std::strcmp(expr, "DENDRO_0931") == 0) {
        if (dir == 0) return grad_0_ccz4_D0931[pp];
        if (dir == 1) return grad_1_ccz4_D0931[pp];
        if (dir == 2) return grad_2_ccz4_D0931[pp];
    }
    if (std::strcmp(expr, "-DENDRO_0293") == 0 && dir == 0) {
        return -grad_0_ccz4_D0293[pp];
    }
    if (std::strcmp(expr, "DENDRO_0521*DENDRO_1021") == 0 && dir == 1) {
        return grad_1_ccz4_D0521_norm[pp];
    }
    if (std::strcmp(expr, "-DENDRO_0293*DENDRO_1021") == 0 && dir == 0) {
        return grad_0_ccz4_mD0293_norm[pp];
    }
    if (std::strcmp(expr, "-DENDRO_0930*DENDRO_1021") == 0 && dir == 2) {
        return grad_2_ccz4_mD0930_norm[pp];
    }
    if (std::strcmp(expr, "DENDRO_0585*DENDRO_1021") == 0 && dir == 0) {
        return grad_0_ccz4_D0585_norm[pp];
    }
    if (std::strcmp(expr, "DENDRO_0426*DENDRO_1021") == 0 && dir == 2) {
        return grad_2_ccz4_D0426_norm[pp];
    }

    std::cerr << "Unlowered CCZ4 grad(" << dir << ", " << expr
              << ") at pp=" << pp << std::endl;
    std::abort();
};

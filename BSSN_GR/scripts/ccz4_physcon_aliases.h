const double p_expo = 1.0;

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
double *const psi4_imag  = psi4_img;

#define CCZ4_PHYSCON_ALIAS_DERIVS(name, base)                                  \
    double *const grad_0_##name    = grad_0_##base;                             \
    double *const grad_1_##name    = grad_1_##base;                             \
    double *const grad_2_##name    = grad_2_##base;                             \
    double *const grad2_0_0_##name = grad2_0_0_##base;                          \
    double *const grad2_0_1_##name = grad2_0_1_##base;                          \
    double *const grad2_0_2_##name = grad2_0_2_##base;                          \
    double *const grad2_1_1_##name = grad2_1_1_##base;                          \
    double *const grad2_1_2_##name = grad2_1_2_##base;                          \
    double *const grad2_2_2_##name = grad2_2_2_##base

CCZ4_PHYSCON_ALIAS_DERIVS(gammat00, gt0);
CCZ4_PHYSCON_ALIAS_DERIVS(gammat01, gt1);
CCZ4_PHYSCON_ALIAS_DERIVS(gammat02, gt2);
CCZ4_PHYSCON_ALIAS_DERIVS(gammat11, gt3);
CCZ4_PHYSCON_ALIAS_DERIVS(gammat12, gt4);
CCZ4_PHYSCON_ALIAS_DERIVS(gammat22, gt5);
#undef CCZ4_PHYSCON_ALIAS_DERIVS

#define CCZ4_PHYSCON_ALIAS_FIRST_DERIVS(name, base)                            \
    double *const grad_0_##name = grad_0_##base;                                \
    double *const grad_1_##name = grad_1_##base;                                \
    double *const grad_2_##name = grad_2_##base

CCZ4_PHYSCON_ALIAS_FIRST_DERIVS(At00, At0);
CCZ4_PHYSCON_ALIAS_FIRST_DERIVS(At01, At1);
CCZ4_PHYSCON_ALIAS_FIRST_DERIVS(At02, At2);
CCZ4_PHYSCON_ALIAS_FIRST_DERIVS(At11, At3);
CCZ4_PHYSCON_ALIAS_FIRST_DERIVS(At12, At4);
CCZ4_PHYSCON_ALIAS_FIRST_DERIVS(At22, At5);
#undef CCZ4_PHYSCON_ALIAS_FIRST_DERIVS

double *const physcon_workspace  = new double[12 * n];
double *const physcon_D0130      = physcon_workspace + 0 * n;
double *const physcon_D0147      = physcon_workspace + 1 * n;
double *const physcon_D0165      = physcon_workspace + 2 * n;
double *const grad_0_physcon_D0130 = physcon_workspace + 3 * n;
double *const grad_1_physcon_D0130 = physcon_workspace + 4 * n;
double *const grad_2_physcon_D0130 = physcon_workspace + 5 * n;
double *const grad_0_physcon_D0147 = physcon_workspace + 6 * n;
double *const grad_1_physcon_D0147 = physcon_workspace + 7 * n;
double *const grad_2_physcon_D0147 = physcon_workspace + 8 * n;
double *const grad_0_physcon_D0165 = physcon_workspace + 9 * n;
double *const grad_1_physcon_D0165 = physcon_workspace + 10 * n;
double *const grad_2_physcon_D0165 = physcon_workspace + 11 * n;

for (unsigned int qq = 0; qq < n; qq++) {
    const unsigned int ii = qq % nx;
    const unsigned int jj = (qq / nx) % ny;
    const unsigned int kk = qq / (nx * ny);
    const double x        = pmin[0] + ii * hx;
    const double y        = pmin[1] + jj * hy;
    const double z        = pmin[2] + kk * hz;
    const unsigned int pp = qq;

#define grad(dir, expr) 0.0
    double DENDRO_0000 = pow(chi[pp], p_expo);
    double DENDRO_0001 = 1.0 / (DENDRO_0000);
    double DENDRO_0002 = x * x;
    double DENDRO_0003 = DENDRO_0002 * gammat00[pp];
    double DENDRO_0004 = y * y;
    double DENDRO_0005 = DENDRO_0004 * gammat11[pp];
    double DENDRO_0006 = z * z;
    double DENDRO_0007 = x * y;
    double DENDRO_0008 = 2 * gammat01[pp];
    double DENDRO_0009 = DENDRO_0007 * DENDRO_0008;
    double DENDRO_0010 = gammat02[pp] * x;
    double DENDRO_0011 = 2 * z;
    double DENDRO_0012 = gammat12[pp] * y;
    double DENDRO_0013 =
        DENDRO_0003 + DENDRO_0005 + DENDRO_0006 * gammat22[pp] +
        DENDRO_0009 + DENDRO_0010 * DENDRO_0011 +
        DENDRO_0011 * DENDRO_0012;
    double DENDRO_0014 = 1.0 / (DENDRO_0013);
    double DENDRO_0015 = DENDRO_0004 * gammat01[pp];
    double DENDRO_0016 = DENDRO_0007 * gammat00[pp];
    double DENDRO_0017 =
        -DENDRO_0002 * gammat01[pp] + DENDRO_0015 + DENDRO_0016 +
        gammat02[pp] * y * z - gammat11[pp] * x * y -
        gammat12[pp] * x * z;
    double DENDRO_0018 = -DENDRO_0017;
    double DENDRO_0019 = DENDRO_0014 * DENDRO_0018;
    double DENDRO_0020 = DENDRO_0002 + DENDRO_0004;
    double DENDRO_0021 = DENDRO_0010 * DENDRO_0020;
    double DENDRO_0022 = DENDRO_0012 * DENDRO_0020;
    double DENDRO_0023 =
        DENDRO_0003 * z + DENDRO_0005 * z +
        DENDRO_0006 * DENDRO_0010 + DENDRO_0006 * DENDRO_0012 +
        DENDRO_0009 * z - DENDRO_0020 * gammat22[pp] * z -
        DENDRO_0021 - DENDRO_0022;
    double DENDRO_0024 = DENDRO_0014 * DENDRO_0023;
    double DENDRO_0025 = -DENDRO_0024 + z;
    double DENDRO_0026 = DENDRO_0020 + DENDRO_0024 * z;
    double DENDRO_0027 =
        -DENDRO_0002 * DENDRO_0025 * gammat01[pp] +
        DENDRO_0015 * DENDRO_0025 + DENDRO_0016 * DENDRO_0025 -
        DENDRO_0025 * gammat11[pp] * x * y -
        DENDRO_0026 * gammat02[pp] * y +
        DENDRO_0026 * gammat12[pp] * x;
    double DENDRO_0028 = -DENDRO_0027;
    double DENDRO_0029 =
        DENDRO_0003 * DENDRO_0006 + DENDRO_0005 * DENDRO_0006 +
        DENDRO_0006 * DENDRO_0009 - DENDRO_0011 * DENDRO_0021 -
        DENDRO_0011 * DENDRO_0022 -
        DENDRO_0014 * DENDRO_0023 * DENDRO_0023 +
        DENDRO_0020 * DENDRO_0020 * gammat22[pp];
    double DENDRO_0030 = 1.0 / (DENDRO_0029);
    double DENDRO_0031 = DENDRO_0001 * DENDRO_0030;
    double DENDRO_0032 =
        DENDRO_0001 * DENDRO_0019 * z -
        DENDRO_0026 * DENDRO_0028 * DENDRO_0031;
    double DENDRO_0033 =
        -DENDRO_0002 * gammat11[pp] - DENDRO_0004 * gammat00[pp] +
        DENDRO_0009;
    double DENDRO_0034 = DENDRO_0014 * DENDRO_0018 * DENDRO_0018 +
                         DENDRO_0028 * DENDRO_0028 * DENDRO_0030 +
                         DENDRO_0033;
    double DENDRO_0035 = 1.0 / (DENDRO_0034);
    double DENDRO_0036 = grad2_2_2_chi[pp];
    double DENDRO_0037 = -DENDRO_0036;
    double DENDRO_0038 = gammat02[pp] * gammat02[pp];
    double DENDRO_0039 = gammat00[pp] * gammat22[pp];
    double DENDRO_0040 = -DENDRO_0038 + DENDRO_0039;
    double DENDRO_0041 = grad_1_gammat22[pp];
    double DENDRO_0042 = 0.5 * DENDRO_0041;
    double DENDRO_0043 = grad_2_gammat12[pp];
    double DENDRO_0044 = 1.0 * DENDRO_0043;
    double DENDRO_0045 = -DENDRO_0044;
    double DENDRO_0046 = DENDRO_0042 + DENDRO_0045;
    double DENDRO_0047 = grad_0_gammat22[pp];
    double DENDRO_0048 = 0.5 * DENDRO_0047;
    double DENDRO_0049 = grad_2_gammat02[pp];
    double DENDRO_0050 = 1.0 * DENDRO_0049;
    double DENDRO_0051 = -DENDRO_0050;
    double DENDRO_0052 = DENDRO_0048 + DENDRO_0051;
    double DENDRO_0053 = gammat02[pp] * gammat12[pp];
    double DENDRO_0054 = -DENDRO_0053 + gammat01[pp] * gammat22[pp];
    double DENDRO_0055 = gammat00[pp] * gammat12[pp] -
                         gammat01[pp] * gammat02[pp];
    double DENDRO_0056 = grad_2_gammat22[pp];
    double DENDRO_0057 = 0.5 * DENDRO_0056;
    double DENDRO_0058 = DENDRO_0055 * DENDRO_0057;
    double DENDRO_0059 =
        DENDRO_0040 * DENDRO_0046 - DENDRO_0052 * DENDRO_0054 +
        DENDRO_0058;
    double DENDRO_0060 = grad_1_chi[pp];
    double DENDRO_0061 = gammat12[pp] * gammat12[pp];
    double DENDRO_0062 = gammat01[pp] * gammat01[pp];
    double DENDRO_0063 =
        -DENDRO_0008 * DENDRO_0053 + DENDRO_0038 * gammat11[pp] -
        DENDRO_0039 * gammat11[pp] + DENDRO_0061 * gammat00[pp] +
        DENDRO_0062 * gammat22[pp];
    double DENDRO_0064 = 1.0 / (DENDRO_0063);
    double DENDRO_0065 = DENDRO_0060 * DENDRO_0064;
    double DENDRO_0066 = -DENDRO_0061 + gammat11[pp] * gammat22[pp];
    double DENDRO_0067 = gammat01[pp] * gammat12[pp] -
                         gammat02[pp] * gammat11[pp];
    double DENDRO_0068 =
        DENDRO_0046 * DENDRO_0054 - DENDRO_0052 * DENDRO_0066 +
        DENDRO_0057 * DENDRO_0067;
    double DENDRO_0069 = -DENDRO_0068;
    double DENDRO_0070 = grad_0_chi[pp];
    double DENDRO_0071 = DENDRO_0064 * DENDRO_0070;
    double DENDRO_0072 = -DENDRO_0062 + gammat00[pp] * gammat11[pp];
    double DENDRO_0073 =
        DENDRO_0046 * DENDRO_0055 - DENDRO_0052 * DENDRO_0067 +
        DENDRO_0057 * DENDRO_0072;
    double DENDRO_0074 = -DENDRO_0073;
    double DENDRO_0075 = grad_2_chi[pp];
    double DENDRO_0076 = DENDRO_0064 * DENDRO_0075;
    double DENDRO_0077 = 1.0 / (chi[pp]);
    double DENDRO_0078 = 6 * DENDRO_0077;
    double DENDRO_0079 = 1.0 / (chi[pp] * chi[pp]);
    double DENDRO_0080 = DENDRO_0075 * DENDRO_0075;
    double DENDRO_0081 = DENDRO_0079 * DENDRO_0080;
    double DENDRO_0082 = 1.0 / (DENDRO_0063 * DENDRO_0063);
    double DENDRO_0083 = DENDRO_0069 * DENDRO_0072;
    double DENDRO_0084 = grad_0_gammat11[pp];
    double DENDRO_0085 = 0.5 * DENDRO_0084;
    double DENDRO_0086 = grad_1_gammat01[pp];
    double DENDRO_0087 = 1.0 * DENDRO_0086;
    double DENDRO_0088 = -DENDRO_0087;
    double DENDRO_0089 = DENDRO_0085 + DENDRO_0088;
    double DENDRO_0090 = grad_1_gammat12[pp];
    double DENDRO_0091 = 1.0 * DENDRO_0090;
    double DENDRO_0092 = grad_2_gammat11[pp];
    double DENDRO_0093 = 0.5 * DENDRO_0092;
    double DENDRO_0094 = DENDRO_0091 - DENDRO_0093;
    double DENDRO_0095 = grad_1_gammat11[pp];
    double DENDRO_0096 = 0.5 * DENDRO_0095;
    double DENDRO_0097 = DENDRO_0054 * DENDRO_0096;
    double DENDRO_0098 =
        DENDRO_0066 * DENDRO_0089 - DENDRO_0067 * DENDRO_0094 +
        DENDRO_0097;
    double DENDRO_0099 = DENDRO_0040 * DENDRO_0098;
    double DENDRO_0100 = grad_0_gammat02[pp];
    double DENDRO_0101 = 1.0 * DENDRO_0100;
    double DENDRO_0102 = grad_2_gammat00[pp];
    double DENDRO_0103 = 0.5 * DENDRO_0102;
    double DENDRO_0104 = DENDRO_0101 - DENDRO_0103;
    double DENDRO_0105 = grad_0_gammat01[pp];
    double DENDRO_0106 = 1.0 * DENDRO_0105;
    double DENDRO_0107 = grad_1_gammat00[pp];
    double DENDRO_0108 = 0.5 * DENDRO_0107;
    double DENDRO_0109 = DENDRO_0106 - DENDRO_0108;
    double DENDRO_0110 = grad_0_gammat00[pp];
    double DENDRO_0111 = 0.5 * DENDRO_0110;
    double DENDRO_0112 =
        -DENDRO_0054 * DENDRO_0109 + DENDRO_0066 * DENDRO_0111 +
        DENDRO_0067 * DENDRO_0104;
    double DENDRO_0113 = -DENDRO_0112;
    double DENDRO_0114 = DENDRO_0066 * DENDRO_0113;
    double DENDRO_0115 = grad_1_gammat02[pp];
    double DENDRO_0116 = grad_0_gammat12[pp];
    double DENDRO_0117 = grad_2_gammat01[pp];
    double DENDRO_0118 = -DENDRO_0115 + DENDRO_0116 + DENDRO_0117;
    double DENDRO_0119 =
        DENDRO_0047 * DENDRO_0067 - DENDRO_0054 * DENDRO_0118 +
        DENDRO_0066 * DENDRO_0102;
    double DENDRO_0120 = -DENDRO_0119;
    double DENDRO_0121 = 1.0 * DENDRO_0067;
    double DENDRO_0122 = DENDRO_0115 - DENDRO_0116 + DENDRO_0117;
    double DENDRO_0123 =
        DENDRO_0041 * DENDRO_0067 - DENDRO_0054 * DENDRO_0092 +
        DENDRO_0066 * DENDRO_0122;
    double DENDRO_0124 = -DENDRO_0123;
    double DENDRO_0125 = DENDRO_0115 + DENDRO_0116 - DENDRO_0117;
    double DENDRO_0126 =
        -DENDRO_0054 * DENDRO_0084 + DENDRO_0066 * DENDRO_0107 +
        DENDRO_0067 * DENDRO_0125;
    double DENDRO_0127 = -DENDRO_0126;
    double DENDRO_0128 =
        -1.0 * DENDRO_0054 * DENDRO_0127 -
        1.0 * DENDRO_0055 * DENDRO_0124 + DENDRO_0083 + DENDRO_0099 +
        DENDRO_0114 + DENDRO_0120 * DENDRO_0121;
    double DENDRO_0129 = -DENDRO_0128;
    double DENDRO_0130 = DENDRO_0082 * DENDRO_0129;
    double DENDRO_0131 = grad(2, DENDRO_0130);
    double DENDRO_0132 = DENDRO_0059 * DENDRO_0072;
    double DENDRO_0133 =
        DENDRO_0040 * DENDRO_0096 + DENDRO_0054 * DENDRO_0089 -
        DENDRO_0055 * DENDRO_0094;
    double DENDRO_0134 = -DENDRO_0133;
    double DENDRO_0135 = DENDRO_0040 * DENDRO_0134;
    double DENDRO_0136 = DENDRO_0054 * DENDRO_0111 +
                         DENDRO_0055 * DENDRO_0104;
    double DENDRO_0137 = -DENDRO_0040 * DENDRO_0109 + DENDRO_0136;
    double DENDRO_0138 = DENDRO_0066 * DENDRO_0137;
    double DENDRO_0139 = DENDRO_0047 * DENDRO_0055 +
                         DENDRO_0054 * DENDRO_0102;
    double DENDRO_0140 = -DENDRO_0040 * DENDRO_0118 + DENDRO_0139;
    double DENDRO_0141 = DENDRO_0041 * DENDRO_0055 +
                         DENDRO_0054 * DENDRO_0122;
    double DENDRO_0142 = -DENDRO_0040 * DENDRO_0092 + DENDRO_0141;
    double DENDRO_0143 = DENDRO_0054 * DENDRO_0107 +
                         DENDRO_0055 * DENDRO_0125;
    double DENDRO_0144 = -DENDRO_0040 * DENDRO_0084 + DENDRO_0143;
    double DENDRO_0145 =
        -1.0 * DENDRO_0054 * DENDRO_0144 -
        1.0 * DENDRO_0055 * DENDRO_0142 + DENDRO_0121 * DENDRO_0140 +
        DENDRO_0132 + DENDRO_0135 + DENDRO_0138;
    double DENDRO_0146 = -DENDRO_0145;
    double DENDRO_0147 = DENDRO_0082 * DENDRO_0146;
    double DENDRO_0148 = grad(2, DENDRO_0147);
    double DENDRO_0149 = 12.0 * gammat12[pp];
    double DENDRO_0150 = DENDRO_0072 * DENDRO_0074;
    double DENDRO_0151 = DENDRO_0055 * DENDRO_0096;
    double DENDRO_0152 =
        DENDRO_0067 * DENDRO_0089 - DENDRO_0072 * DENDRO_0094 +
        DENDRO_0151;
    double DENDRO_0153 = DENDRO_0040 * DENDRO_0152;
    double DENDRO_0154 =
        -DENDRO_0055 * DENDRO_0109 + DENDRO_0067 * DENDRO_0111 +
        DENDRO_0072 * DENDRO_0104;
    double DENDRO_0155 = -DENDRO_0154;
    double DENDRO_0156 = DENDRO_0066 * DENDRO_0155;
    double DENDRO_0157 =
        DENDRO_0047 * DENDRO_0072 - DENDRO_0055 * DENDRO_0118 +
        DENDRO_0067 * DENDRO_0102;
    double DENDRO_0158 = -DENDRO_0157;
    double DENDRO_0159 =
        DENDRO_0041 * DENDRO_0072 - DENDRO_0055 * DENDRO_0092 +
        DENDRO_0067 * DENDRO_0122;
    double DENDRO_0160 = -DENDRO_0159;
    double DENDRO_0161 =
        -DENDRO_0055 * DENDRO_0084 + DENDRO_0067 * DENDRO_0107 +
        DENDRO_0072 * DENDRO_0125;
    double DENDRO_0162 = -DENDRO_0161;
    double DENDRO_0163 =
        -1.0 * DENDRO_0054 * DENDRO_0162 -
        1.0 * DENDRO_0055 * DENDRO_0160 + DENDRO_0121 * DENDRO_0158 +
        DENDRO_0150 + DENDRO_0153 + DENDRO_0156;
    double DENDRO_0164 = -DENDRO_0163;
    double DENDRO_0165 = DENDRO_0082 * DENDRO_0164;
#undef grad

    physcon_D0130[qq] = DENDRO_0130;
    physcon_D0147[qq] = DENDRO_0147;
    physcon_D0165[qq] = DENDRO_0165;
}

deriv_x(grad_0_physcon_D0130, physcon_D0130, hx, sz, bflag);
deriv_y(grad_1_physcon_D0130, physcon_D0130, hy, sz, bflag);
deriv_z(grad_2_physcon_D0130, physcon_D0130, hz, sz, bflag);
deriv_x(grad_0_physcon_D0147, physcon_D0147, hx, sz, bflag);
deriv_y(grad_1_physcon_D0147, physcon_D0147, hy, sz, bflag);
deriv_z(grad_2_physcon_D0147, physcon_D0147, hz, sz, bflag);
deriv_x(grad_0_physcon_D0165, physcon_D0165, hx, sz, bflag);
deriv_y(grad_1_physcon_D0165, physcon_D0165, hy, sz, bflag);
deriv_z(grad_2_physcon_D0165, physcon_D0165, hz, sz, bflag);

auto physcon_grad = [&](const int dir, const char *expr,
                        const unsigned int pp) -> double {
    if (std::strcmp(expr, "DENDRO_0130") == 0) {
        if (dir == 0) return grad_0_physcon_D0130[pp];
        if (dir == 1) return grad_1_physcon_D0130[pp];
        if (dir == 2) return grad_2_physcon_D0130[pp];
    }
    if (std::strcmp(expr, "DENDRO_0147") == 0) {
        if (dir == 0) return grad_0_physcon_D0147[pp];
        if (dir == 1) return grad_1_physcon_D0147[pp];
        if (dir == 2) return grad_2_physcon_D0147[pp];
    }
    if (std::strcmp(expr, "DENDRO_0165") == 0) {
        if (dir == 0) return grad_0_physcon_D0165[pp];
        if (dir == 1) return grad_1_physcon_D0165[pp];
        if (dir == 2) return grad_2_physcon_D0165[pp];
    }
    std::cerr << "Unlowered CCZ4 physcon grad(" << dir << ", " << expr
              << ") at pp=" << pp << std::endl;
    std::abort();
};

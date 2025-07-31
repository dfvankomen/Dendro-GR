"""
BSSN core variables .
"""

import argparse
import sys as sys

import dendro
from sympy import *

###################################################################
# initialize
###################################################################

l1, l2, l3, l4, eta = symbols("lambda[0] lambda[1] lambda[2] lambda[3] eta")
lf0, lf1 = symbols("lambda_f[0] lambda_f[1]")
kappa = symbols("kappa")
# Additional parameters for damping term
R0 = symbols("BSSN_ETA_R0")
ep1, ep2 = symbols("BSSN_ETA_POWER[0] BSSN_ETA_POWER[1]")

xi1, xi2, xi3 = symbols("BSSN_XI[0] BSSN_XI[1] BSSN_XI[2] ")

# declare variables
a = dendro.scalar("alpha", "[pp]")
chi = dendro.scalar("chi", "[pp]")
K = dendro.scalar("K", "[pp]")

Gt = dendro.vec3("Gt", "[pp]")
b = dendro.vec3("beta", "[pp]")
B = dendro.vec3("B", "[pp]")

gt = dendro.sym_3x3("gt", "[pp]")
At = dendro.sym_3x3("At", "[pp]")

Gt_rhs = dendro.vec3("Gt_rhs", "[pp]")

PerpE = dendro.vec3("PerpE", "[pp]")

PerpB = dendro.vec3("PerpB", "[pp]")

DampPsi = dendro.scalar("DampPsi", "[pp]")

DampPhi =dendro.scalar("DampPhi", "[pp]")

# Lie derivative weight
weight = -Rational(2, 3)
weight_Gt = Rational(2, 3)

# standard Levi-Civita symbol
leviCivita = [
    [[0, 0, 0], [0, 0, 1], [0, -1, 0]],
    [[0, 0, -1], [0, 0, 0], [1, 0, 0]],
    [[0, 1, 0], [-1, 0, 0], [0, 0, 0]],
]

# specify the functions for computing first and second derivatives

# d = dendro.set_first_derivative('grad')    # first argument is direction
# d2s = dendro.set_second_derivative('grad2')  # first 2 arguments are directions
# ad = dendro.set_advective_derivative('agrad')  # first argument is direction
# kod = dendro.set_kreiss_oliger_dissipation('kograd')
dendro.d = lambda i, x: symbols("grad_%d_%s" % (i, x))
dendro.d2 = lambda i, j, x: symbols("grad2_%d_%d_%s" % (min(i, j), max(i, j), x))
dendro.ad = dendro.d
dendro.kod = dendro.undef

d = dendro.d
ad = dendro.ad
kod = dendro.kod
d2 = dendro.d2

t = symbols("t")  # time; needed for SSL

# add symbols used for CAHD
ham = symbols("ham[pp]")  # hamiltonian constraint violation
C_CAHD = symbols("BSSN_CAHD_C")  # coefficient for CAHD strength
dt = symbols("dt")  # simulation time step
dx_i = symbols("dx_i")  # spatial resolution of current grid
dx_min = symbols("dx_min")  # spatial resolution of finest grid

dendro.set_metric(gt)
igt = dendro.get_inverse_metric()

ham_temp_var = symbols("ham_temp")


eta_func = (
    R0
    * sqrt(sum([igt[i, j] * d(i, chi) * d(j, chi) for i, j in dendro.e_ij]))
    / ((1 - chi**ep1) ** ep2)
)


def bssn_puncture_gauge(
    eta_damp, isStaged=False, prefix="", sslGaugeCondition=False, enableCAHD=False
):
        # set up all 10 parts of the stress energy tensor



    perpTpart4 = (
        1
        / (8 * pi)
        * Matrix(
            [
                2
                / (pow(chi, 2.0))
                * sum(
                    [
                        PerpB[k]
                        * PerpB[m]
                        * leviCivita[i][j][k]
                        * leviCivita[l][f][m]
                        * igt[j, f]
                        for k in dendro.e_i
                        for j in dendro.e_i
                        for m in dendro.e_i
                        for f in dendro.e_i
                    ]
                )
                for i, l in dendro.e_ij
            ]
        )
    )

    perpTpart5 = (
        1
        / (8 * pi)
        * Matrix(
            [
                -2
                / (chi**2)
                * sum(
                    [
                        gt[i, k] * gt[j, l] * PerpE[k] * PerpE[l]
                        for k in dendro.e_i
                        for l in dendro.e_i
                    ]
                )
                for i, j in dendro.e_ij
            ]
        )
    )
    perpTpart6 = (
        1
        / (8 * pi)
        * Matrix(
            [
                1
                / (chi**2)
                * sum(
                    [
                        PerpB[l] * PerpB[m] * gt[l, m]
                        for l in dendro.e_i
                        for m in dendro.e_i
                    ]
                )
                * gt[i, j]
                for i, j in dendro.e_ij
            ]
        )
    )
    perpTpart7 = (
        1
        / (8 * pi)
        * Matrix(
            [

                1/ (chi**2)
                * sum(
                    [
                        PerpE[l] * PerpE[m] * gt[l, m]
                        for l in dendro.e_i
                        for m in dendro.e_i
                    ]
                )
                * gt[i, j]
                for i, j in dendro.e_ij
            ]
        )
    )
   
   
    perpTpart4 = perpTpart4.reshape(3, 3)
    perpTpart5 = perpTpart5.reshape(3, 3)
    perpTpart6 = perpTpart6.reshape(3, 3)
    perpTpart7 = perpTpart7.reshape(3, 3)
   

    # assemble stress energy tensor, shows up in At_rhs
    perpT = Matrix(
        [

            + perpTpart4[i, j]
            + perpTpart5[i, j]
            + perpTpart6[i, j]
            + perpTpart7[i, j]

            for i, j in dendro.e_ij
        ]
    )
    perpT = perpT.reshape(
        3, 3
    )  # very key. Before, it would return a 9x1 sym.Matrix I think
    # stress energy density, shows up in trK_rhs
    rho = (
        1
        / (8 * pi)
        * (sum(
                [
                    PerpE[i] * PerpE[j] * gt[i, j]
                    for i in dendro.e_i
                    for j in dendro.e_i
                ]
            )
            / chi
            + 
            sum(
                [
                    PerpB[i] * PerpB[j] * gt[i, j]
                    for i in dendro.e_i
                    for j in dendro.e_i
                ]
            )
            / chi
        )
    )
    # stress energy current, shows up in Gt_rhs
    stressCurrent = [
        1
        / (8 * pi)
        * 2
        * (chi) ** (Rational(1, 2))
        * sum(
            [
                leviCivita[i][j][k] * gt[j, l] * PerpE[l] * gt[k, m] * PerpB[m]
                for j in dendro.e_i
                for k in dendro.e_i
                for l in dendro.e_i
                for m in dendro.e_i
            ]
        )
        for i in dendro.e_i
    ]
# compute the trace of the stress energy tensor, this shows up in trK_rhs
    perpTTrace = dendro.trace(perpT)
# [EWH]:  I don't think these two quantities (the covariant version of the
#         normal vector to the spacelike hypersurfaces and the mean extrinsic
#         curvature used in the gauge conditions) get used, so I will comment
#         them out.
# nc = (-a, 0, 0, 0)
# trK0 = 0
    if not isStaged:
        C1 = dendro.get_first_christoffel()
        C2 = dendro.get_second_christoffel()
        C2_spatial = dendro.get_complete_christoffel(chi)
        [R, Rt, Rphi, CalGt] = dendro.compute_ricci(Gt, chi)

        ham_computation = (
            sum(chi * igt[j, k] * R[j, k] for j, k in dendro.e_ij)
            - dendro.sqr(At)
            + Rational(2, 3) * K**2
        )

        if sslGaugeCondition:
            # enable slow-start lapse
            W = chi**0.5

            h = symbols("h_ssl")
            sig = symbols("sig_ssl")
            # h = 0.6
            # sig = 20
            a_rhs = (
                l1 * dendro.lie(b, a)
                - 2 * a * K
                - W * (h * exp(-(t**2) / (2 * sig**2))) * (a - W)
            )
        else:
            a_rhs = l1 * dendro.lie(b, a) - 2 * a * K

        b_rhs = [
            (Rational(3, 4) * (lf0 + lf1 * a) * B[i] + l2 * dendro.vec_j_ad_j(b, b[i]))
            for i in dendro.e_i
        ]

        gt_rhs = dendro.lie(b, gt, weight) - 2 * a * At

        chi_rhs = dendro.lie(b, chi, weight) + Rational(2, 3) * (chi * a * K)

        if enableCAHD:
            # turn on curvature-adjusted Hamiltonian-constraint damping
            # chi_rhs += C_CAHD * chi * (dt * dx_i / dx_min) * ham # Etienne's method
            chi_rhs += C_CAHD * chi * (dx_i**2 / dt) * ham_computation  # WKB's method

        AikAkj = Matrix(
            [
                sum(
                    [
                        At[ii, kk]
                        * sum([igt[kk, ll] * At[ll, jj] for ll in dendro.e_i])
                        for kk in dendro.e_i
                    ]
                )
                for ii, jj in dendro.e_ij
            ]
        ).reshape(3, 3)

        At_rhs = (
            (dendro.lie(b, At, weight))
            + chi
            * dendro.trace_free(
                a * R - dendro.DiDj(a) - 8 * pi * a * perpT
            )
            + a * (K * At - 2 * AikAkj)
        )

        K_rhs = (
            dendro.lie(b, K)
            - dendro.laplacian(a, chi)
            + a * (K * K / 3 + dendro.sqr(At))
            + 4 * pi * a * (rho + perpTTrace)
        )

        At_UU = dendro.up_up(At)

        Gt_rhs = (
            Matrix([sum(b[j] * ad(j, Gt[i]) for j in dendro.e_i) for i in dendro.e_i])
            - Matrix(
                [sum(CalGt[j] * d(j, b[i]) for j in dendro.e_i) for i in dendro.e_i]
            )
            + Rational(2, 3)
            * Matrix(
                [CalGt[i] * sum(d(j, b[j]) for j in dendro.e_i) for i in dendro.e_i]
            )
            + Matrix(
                [
                    sum(
                        [
                            igt[j, k] * d2(j, k, b[i]) + igt[i, j] * d2(j, k, b[k]) / 3
                            for j, k in dendro.e_ij
                        ]
                    )
                    for i in dendro.e_i
                ]
            )
            - Matrix(
                [
                    sum([2 * At_UU[i, j] * d(j, a) for j in dendro.e_i])
                    for i in dendro.e_i
                ]
            )
            + Matrix(
                [
                    sum(
                        [
                            2 * a * dendro.C2[i, j, k] * At_UU[j, k]
                            for j, k in dendro.e_ij
                        ]
                    )
                    for i in dendro.e_i
                ]
            )
            - Matrix(
                [
                    sum(
                        [
                            a
                            * (
                                3 / chi * At_UU[i, j] * d(j, chi)
                                + Rational(4, 3) * dendro.inv_metric[i, j] * d(j, K)
                            )
                            for j in dendro.e_i
                        ]
                    )
                    for i in dendro.e_i
                ]
            )
            
                +Matrix(
            ([-16 * pi * a/ chi * stressCurrent[i] for i in dendro.e_i])
        )
        ) 

        Gt_rhs = [item for sublist in Gt_rhs.tolist() for item in sublist]

        B_rhs = [
            (
                Gt_rhs[i]
                - eta_damp * B[i]
                + l3 * dendro.vec_j_ad_j(b, B[i])
                - l4 * dendro.vec_j_ad_j(b, Gt[i])
            )
            for i in dendro.e_i
        ]
        PerpE_rhs = [
            dendro.lie(b, PerpE)[i]
            + a * K * PerpE[i]
            - a
            * chi
            * sum([igt[i, j] * d(j, DampPsi) for j in dendro.e_i])
            - pow(chi, Rational(1, 2))
            * sum(
                [
                    leviCivita[j][i][k] * gt[k, l] * PerpB[l] * d(j, a)
                    for j in dendro.e_i
                    for k in dendro.e_i
                    for l in dendro.e_i
                ]
            )
            - a
            * pow(chi, Rational(1, 2))
            * sum(
                [
                    leviCivita[j][i][k] * PerpB[l] * d(j, gt[k, l])
                    for j in dendro.e_i
                    for k in dendro.e_i
                    for l in dendro.e_i
                ]
            )
            - a
            * pow(chi, Rational(1, 2))
            * sum(
                [
                    leviCivita[j][i][k] * gt[k, l] * d(j, PerpB[l])
                    for j in dendro.e_i
                    for k in dendro.e_i
                    for l in dendro.e_i
                ]
            )
            + a
            * pow(chi, -Rational(1, 2))
            * sum(
                [
                    leviCivita[j][i][k] * gt[k, l] * PerpB[l] * d(j, chi)
                    for j in dendro.e_i
                    for k in dendro.e_i
                    for l in dendro.e_i
                ]
            )
            for i in dendro.e_i
        ]


        PerpB_rhs = [
            a * K * PerpB[i]
            - a
            * chi
            * sum([igt[i, j] * d(j, DampPhi) for j in dendro.e_i])
            - pow(chi, Rational(1, 2))
            * sum(
                [
                    leviCivita[i][j][k] * gt[k, l] * PerpE[l] * d(j, a)
                    for j in dendro.e_i
                    for k in dendro.e_i
                    for l in dendro.e_i
                ]
            )
            + pow(chi, -Rational(1, 2))
            * sum(
                [
                    leviCivita[i][j][k] * a * gt[k, l] * PerpE[l] * d(j, chi)
                    for j in dendro.e_i
                    for k in dendro.e_i
                    for l in dendro.e_i
                ]
            )
            - a
            * pow(chi, Rational(1, 2))
            * sum(
                [
                    leviCivita[i][j][k] * PerpE[l] * d(j, gt[k, l])
                    for j in dendro.e_i
                    for k in dendro.e_i
                    for l in dendro.e_i
                ]
            )
            - a
            * pow(chi, Rational(1, 2))
            * sum(
                [
                    leviCivita[i][j][k] * gt[k, l] * d(j, PerpE[l])
                    for j in dendro.e_i
                    for k in dendro.e_i
                    for l in dendro.e_i
                ]
            )
            + dendro.lie(b, PerpB)[i]
            for i in dendro.e_i
        ]

        DampPhi_rhs = (
                        dendro.vec_j_del_j(b, DampPhi)
            - a * sum([d(i, PerpB[i]) for i in dendro.e_i])
            + 3
            * a
            / (2 * chi)
            * sum([PerpB[i] * d(i, chi) for i in dendro.e_i])
            - a * 0.4 * DampPhi
        )

        DampPsi_rhs = (
            dendro.vec_j_del_j(b, DampPsi)
            - a * sum([d(i, PerpE[i]) for i in dendro.e_i])
            + 3
            * a
            / (2 * chi)
            * sum([PerpE[i] * d(i, chi) for i in dendro.e_i])
            - a * 0.4 * DampPsi
        )

        ###################################################################
        # generate code
        ###################################################################

        outs = [a_rhs, b_rhs, gt_rhs, chi_rhs, At_rhs, K_rhs, Gt_rhs, B_rhs, PerpE_rhs, PerpB_rhs, DampPhi_rhs, DampPsi_rhs ]
        vnames = [
            "a_rhs",
            "b_rhs",
            "gt_rhs",
            "chi_rhs",
            "At_rhs",
            "K_rhs",
            "Gt_rhs",
            "B_rhs",
            "PerpE_rhs",
            "PerpB_rhs",
            "DampPhi_rhs",
            "DampPsi_rhs",
        ]
        dendro.generate_cpu(outs, vnames, "[pp]")

    else:
        # note: these are just the symbolic vars that is being used to generate the
        # Gt_rhs by satges

        _Gt_rhs_s1 = dendro.vec3("Gt_rhs_s1_", "[pp]")
        _Gt_rhs_s2 = dendro.vec3("Gt_rhs_s2_", "[pp]")
        _Gt_rhs_s3 = dendro.vec3("Gt_rhs_s3_", "[pp]")
        _Gt_rhs_s4 = dendro.vec3("Gt_rhs_s4_", "[pp]")
        _Gt_rhs_s5 = dendro.vec3("Gt_rhs_s5_", "[pp]")
        _Gt_rhs_s6 = dendro.vec3("Gt_rhs_s6_", "[pp]")
        _Gt_rhs_s7 = dendro.vec3("Gt_rhs_s7_", "[pp]")
        _CalGt = dendro.vec3("CalGt", "[pp]")
        _Gt_rhs = dendro.vec3("Gt_rhs", "[pp]")

        # Gt_rhs staged vars that is being used to generate the code.
        At_UU = dendro.sym_3x3("At_UU", "[pp]")
        CalGt = dendro.vec3("CalGt", "[pp]")
        Gt_rhs_s1 = dendro.vec3("Gt_rhs_s1_", "[pp]")
        Gt_rhs_s2 = dendro.vec3("Gt_rhs_s2_", "[pp]")
        Gt_rhs_s3 = dendro.vec3("Gt_rhs_s3_", "[pp]")
        Gt_rhs_s4 = dendro.vec3("Gt_rhs_s4_", "[pp]")
        Gt_rhs_s5 = dendro.vec3("Gt_rhs_s5_", "[pp]")
        Gt_rhs_s6 = dendro.vec3("Gt_rhs_s6_", "[pp]")
        Gt_rhs_s7 = dendro.vec3("Gt_rhs_s7_", "[pp]")

        C1 = dendro.get_first_christoffel()
        C2 = dendro.get_second_christoffel()
        C2_spatial = dendro.get_complete_christoffel(chi)
        [R, Rt, Rphi, CalGt] = dendro.compute_ricci(Gt, chi)

        if sslGaugeCondition:
            W = chi**0.5
            h = 0.6
            sig = 20
            a_rhs = (
                l1 * dendro.lie(b, a)
                - 2 * a * K
                - W * (h * exp(-(t**2) / (2 * sig**2))) * (a - W)
            )
        else:
            a_rhs = l1 * dendro.lie(b, a) - 2 * a * K

        b_rhs = [
            (Rational(3, 4) * (lf0 + lf1 * a) * B[i] + l2 * dendro.vec_j_ad_j(b, b[i]))
            for i in dendro.e_i
        ]

        gt_rhs = dendro.lie(b, gt, weight) - 2 * a * At

        chi_rhs = dendro.lie(b, chi, weight) + Rational(2, 3) * (chi * a * K)

        if enableCAHD:
            # turn on curvature-adjusted Hamiltonian-constraint damping
            # chi_rhs += C_CAHD * chi * (dt * dx_i / dx_min) * ham # Etienne's method
            chi_rhs += C_CAHD * chi * (dx_i**2 / dt) * ham  # WKB's method

        AikAkj = Matrix(
            [
                sum(
                    [
                        At[i, k]
                        * sum([dendro.inv_metric[k, l] * At[l, j] for l in dendro.e_i])
                        for k in dendro.e_i
                    ]
                )
                for i, j in dendro.e_ij
            ]
        )

        At_rhs = (
            dendro.lie(b, At, weight)
            + chi * dendro.trace_free(a * R - dendro.DiDj(a))
            + a * (K * At - 2 * AikAkj.reshape(3, 3))
        )

        K_rhs = (
            dendro.lie(b, K)
            - dendro.laplacian(a, chi)
            + a * (K * K / 3 + dendro.sqr(At))
        )

        At_UU = dendro.up_up(At)

        Gt_rhs_s1 = [sum(b[j] * ad(j, Gt[i]) for j in dendro.e_i) for i in dendro.e_i]
        Gt_rhs_s2 = [
            sum(_CalGt[j] * d(j, b[i]) for j in dendro.e_i) for i in dendro.e_i
        ]
        Gt_rhs_s3 = [
            _CalGt[i] * sum(d(j, b[j]) for j in dendro.e_i) for i in dendro.e_i
        ]
        Gt_rhs_s4 = [
            sum(
                [
                    igt[j, k] * d2(j, k, b[i]) + igt[i, j] * d2(j, k, b[k]) / 3
                    for j, k in dendro.e_ij
                ]
            )
            for i in dendro.e_i
        ]
        Gt_rhs_s5 = [
            sum([2 * At_UU[i, j] * d(j, a) for j in dendro.e_i]) for i in dendro.e_i
        ]
        Gt_rhs_s6 = [
            sum([2 * a * dendro.C2[i, j, k] * At_UU[j, k] for j, k in dendro.e_ij])
            for i in dendro.e_i
        ]
        Gt_rhs_s7 = [
            sum(
                [
                    a
                    * (
                        3 / chi * At_UU[i, j] * d(j, chi)
                        + Rational(4, 3) * dendro.inv_metric[i, j] * d(j, K)
                    )
                    for j in dendro.e_i
                ]
            )
            for i in dendro.e_i
        ]

        Gt_rhs = (
            Matrix(_Gt_rhs_s1)
            - Matrix(_Gt_rhs_s2)
            + Rational(2, 3) * Matrix(_Gt_rhs_s3)
            + Matrix(_Gt_rhs_s4)
            - Matrix(_Gt_rhs_s5)
            + Matrix(_Gt_rhs_s6)
            - Matrix(_Gt_rhs_s7)
        )

        Gt_rhs = [item for sublist in Gt_rhs.tolist() for item in sublist]

        B_rhs = [
            (
                Gt_rhs[i]
                - eta_damp * B[i]
                + l3 * dendro.vec_j_ad_j(b, B[i])
                - l4 * dendro.vec_j_ad_j(b, Gt[i])
            )
            for i in dendro.e_i
        ]

        outs = [
            a_rhs,
            b_rhs,
            gt_rhs,
            chi_rhs,
            At_rhs,
            K_rhs,
            CalGt,
            Gt_rhs_s1,
            Gt_rhs_s2,
            Gt_rhs_s3,
            Gt_rhs_s4,
            Gt_rhs_s5,
            Gt_rhs_s6,
            Gt_rhs_s7,
            Gt_rhs,
            B_rhs,
        ]
        vnames = [
            "a_rhs",
            "b_rhs",
            "gt_rhs",
            "chi_rhs",
            "At_rhs",
            "K_rhs",
            "CalGt",
            "Gt_rhs_s1_",
            "Gt_rhs_s2_",
            "Gt_rhs_s3_",
            "Gt_rhs_s4_",
            "Gt_rhs_s5_",
            "Gt_rhs_s6_",
            "Gt_rhs_s7_",
            "Gt_rhs",
            "B_rhs",
        ]

        ###################################################################
        # generate code
        ###################################################################

        numVars = len(outs)
        for i in range(0, numVars):
            dendro.generate_separate([outs[i]], [vnames[i]], "[pp]")


def bssn_rochester_puncture_gauge(
    eta_damp, isStaged=False, prefix="", sslGaugeCondition=False
):
    """
    Uses Rochester puncture gauge.
    """

    if not isStaged:
        C1 = dendro.get_first_christoffel()
        C2 = dendro.get_second_christoffel()
        C2_spatial = dendro.get_complete_christoffel(chi)
        [R, Rt, Rphi, CalGt] = dendro.compute_ricci(Gt, chi)

        if sslGaugeCondition:
            # enable slow-start lapse
            W = chi**0.5
            h = 0.6  # M; Gaussian height
            sig = 20  # M; Gaussian stddev
            a_rhs = (
                l1 * dendro.lie(b, a)
                - 2 * a * K
                - W * (h * exp(-(t**2) / (2 * sig**2))) * (a - W)
            )
        else:  # no SSL
            a_rhs = l1 * dendro.lie(b, a) - 2 * a * K

        b_rhs = [
            (
                xi2 * dendro.vec_j_ad_j(b, b[i])
                + Rational(3, 4) * xi3 * Gt[i]
                - eta_damp * b[i]
            )
            for i in dendro.e_i
        ]

        gt_rhs = dendro.lie(b, gt, weight) - 2 * a * At

        chi_rhs = dendro.lie(b, chi, weight) + Rational(2, 3) * (chi * a * K)

        AikAkj = Matrix(
            [
                sum(
                    [
                        At[i, k]
                        * sum([dendro.inv_metric[k, l] * At[l, j] for l in dendro.e_i])
                        for k in dendro.e_i
                    ]
                )
                for i, j in dendro.e_ij
            ]
        )

        At_rhs = (
            dendro.lie(b, At, weight)
            + chi * dendro.trace_free(a * R - dendro.DiDj(a))
            + a * (K * At - 2 * AikAkj.reshape(3, 3))
        )

        K_rhs = (
            dendro.lie(b, K)
            - dendro.laplacian(a, chi)
            + a * (K * K / 3 + dendro.sqr(At))
        )

        At_UU = dendro.up_up(At)

        B_rhs = 0

        Gt_rhs = (
            Matrix([sum(b[j] * ad(j, Gt[i]) for j in dendro.e_i) for i in dendro.e_i])
            - Matrix(
                [sum(CalGt[j] * d(j, b[i]) for j in dendro.e_i) for i in dendro.e_i]
            )
            + Rational(2, 3)
            * Matrix(
                [CalGt[i] * sum(d(j, b[j]) for j in dendro.e_i) for i in dendro.e_i]
            )
            + Matrix(
                [
                    sum(
                        [
                            igt[j, k] * d2(j, k, b[i]) + igt[i, j] * d2(j, k, b[k]) / 3
                            for j, k in dendro.e_ij
                        ]
                    )
                    for i in dendro.e_i
                ]
            )
            - Matrix(
                [
                    sum([2 * At_UU[i, j] * d(j, a) for j in dendro.e_i])
                    for i in dendro.e_i
                ]
            )
            + Matrix(
                [
                    sum(
                        [
                            2 * a * dendro.C2[i, j, k] * At_UU[j, k]
                            for j, k in dendro.e_ij
                        ]
                    )
                    for i in dendro.e_i
                ]
            )
            - Matrix(
                [
                    sum(
                        [
                            a
                            * (
                                3 / chi * At_UU[i, j] * d(j, chi)
                                + Rational(4, 3) * dendro.inv_metric[i, j] * d(j, K)
                            )
                            for j in dendro.e_i
                        ]
                    )
                    for i in dendro.e_i
                ]
            )
        )

        Gt_rhs = [item for sublist in Gt_rhs.tolist() for item in sublist]

        ###################################################################
        # generate code
        ###################################################################

        outs = [a_rhs, b_rhs, gt_rhs, chi_rhs, At_rhs, K_rhs, Gt_rhs]
        vnames = ["a_rhs", "b_rhs", "gt_rhs", "chi_rhs", "At_rhs", "K_rhs", "Gt_rhs"]
        dendro.generate_cpu(outs, vnames, "[pp]")

    else:
        _Gt_rhs_s1 = dendro.vec3("Gt_rhs_s1_", "[pp]")
        _Gt_rhs_s2 = dendro.vec3("Gt_rhs_s2_", "[pp]")
        _Gt_rhs_s3 = dendro.vec3("Gt_rhs_s3_", "[pp]")
        _Gt_rhs_s4 = dendro.vec3("Gt_rhs_s4_", "[pp]")
        _Gt_rhs_s5 = dendro.vec3("Gt_rhs_s5_", "[pp]")
        _Gt_rhs_s6 = dendro.vec3("Gt_rhs_s6_", "[pp]")
        _Gt_rhs_s7 = dendro.vec3("Gt_rhs_s7_", "[pp]")
        _CalGt = dendro.vec3("CalGt", "[pp]")
        _Gt_rhs = dendro.vec3("Gt_rhs", "[pp]")

        # Gt_rhs staged vars that is being used to generate the code.
        At_UU = dendro.sym_3x3("At_UU", "[pp]")
        CalGt = dendro.vec3("CalGt", "[pp]")
        Gt_rhs_s1 = dendro.vec3("Gt_rhs_s1_", "[pp]")
        Gt_rhs_s2 = dendro.vec3("Gt_rhs_s2_", "[pp]")
        Gt_rhs_s3 = dendro.vec3("Gt_rhs_s3_", "[pp]")
        Gt_rhs_s4 = dendro.vec3("Gt_rhs_s4_", "[pp]")
        Gt_rhs_s5 = dendro.vec3("Gt_rhs_s5_", "[pp]")
        Gt_rhs_s6 = dendro.vec3("Gt_rhs_s6_", "[pp]")
        Gt_rhs_s7 = dendro.vec3("Gt_rhs_s7_", "[pp]")

        C1 = dendro.get_first_christoffel()
        C2 = dendro.get_second_christoffel()
        C2_spatial = dendro.get_complete_christoffel(chi)
        [R, Rt, Rphi, CalGt] = dendro.compute_ricci(Gt, chi)

        if sslGaugeCondition:
            # enable slow-start lapse
            W = chi**0.5
            h = 0.6  # M; Gaussian height
            sig = 20  # M; Gaussian stddev
            a_rhs = (
                l1 * dendro.lie(b, a)
                - 2 * a * K
                - W * (h * exp(-(t**2) / (2 * sig**2))) * (a - W)
            )
        else:  # no SSL
            a_rhs = l1 * dendro.lie(b, a) - 2 * a * K

        b_rhs = [
            (
                xi2 * dendro.vec_j_ad_j(b, b[i])
                + Rational(3, 4) * xi3 * Gt[i]
                - eta_damp * b[i]
            )
            for i in dendro.e_i
        ]

        gt_rhs = dendro.lie(b, gt, weight) - 2 * a * At

        chi_rhs = dendro.lie(b, chi, weight) + Rational(2, 3) * (chi * a * K)

        AikAkj = Matrix(
            [
                sum(
                    [
                        At[i, k]
                        * sum([dendro.inv_metric[k, l] * At[l, j] for l in dendro.e_i])
                        for k in dendro.e_i
                    ]
                )
                for i, j in dendro.e_ij
            ]
        )

        At_rhs = (
            dendro.lie(b, At, weight)
            + chi * dendro.trace_free(a * R - dendro.DiDj(a))
            + a * (K * At - 2 * AikAkj.reshape(3, 3))
        )

        K_rhs = (
            dendro.lie(b, K)
            - dendro.laplacian(a, chi)
            + a * (K * K / 3 + dendro.sqr(At))
        )

        At_UU = dendro.up_up(At)

        B_rhs = 0

        Gt_rhs_s1 = [sum(b[j] * ad(j, Gt[i]) for j in dendro.e_i) for i in dendro.e_i]
        Gt_rhs_s2 = [
            sum(_CalGt[j] * d(j, b[i]) for j in dendro.e_i) for i in dendro.e_i
        ]
        Gt_rhs_s3 = [
            _CalGt[i] * sum(d(j, b[j]) for j in dendro.e_i) for i in dendro.e_i
        ]
        Gt_rhs_s4 = [
            sum(
                [
                    igt[j, k] * d2(j, k, b[i]) + igt[i, j] * d2(j, k, b[k]) / 3
                    for j, k in dendro.e_ij
                ]
            )
            for i in dendro.e_i
        ]
        Gt_rhs_s5 = [
            sum([2 * At_UU[i, j] * d(j, a) for j in dendro.e_i]) for i in dendro.e_i
        ]
        Gt_rhs_s6 = [
            sum([2 * a * dendro.C2[i, j, k] * At_UU[j, k] for j, k in dendro.e_ij])
            for i in dendro.e_i
        ]
        Gt_rhs_s7 = [
            sum(
                [
                    a
                    * (
                        3 / chi * At_UU[i, j] * d(j, chi)
                        + Rational(4, 3) * dendro.inv_metric[i, j] * d(j, K)
                    )
                    for j in dendro.e_i
                ]
            )
            for i in dendro.e_i
        ]

        Gt_rhs = (
            Matrix(_Gt_rhs_s1)
            - Matrix(_Gt_rhs_s2)
            + Rational(2, 3) * Matrix(_Gt_rhs_s3)
            + Matrix(_Gt_rhs_s4)
            - Matrix(_Gt_rhs_s5)
            + Matrix(_Gt_rhs_s6)
            - Matrix(_Gt_rhs_s7)
        )

        Gt_rhs = [item for sublist in Gt_rhs.tolist() for item in sublist]

        outs = [
            a_rhs,
            b_rhs,
            gt_rhs,
            chi_rhs,
            At_rhs,
            K_rhs,
            CalGt,
            Gt_rhs_s1,
            Gt_rhs_s2,
            Gt_rhs_s3,
            Gt_rhs_s4,
            Gt_rhs_s5,
            Gt_rhs_s6,
            Gt_rhs_s7,
            Gt_rhs,
        ]
        vnames = [
            "a_rhs",
            "b_rhs",
            "gt_rhs",
            "chi_rhs",
            "At_rhs",
            "K_rhs",
            "CalGt",
            "Gt_rhs_s1_",
            "Gt_rhs_s2_",
            "Gt_rhs_s3_",
            "Gt_rhs_s4_",
            "Gt_rhs_s5_",
            "Gt_rhs_s6_",
            "Gt_rhs_s7_",
            "Gt_rhs",
        ]

        ###################################################################
        # generate code
        ###################################################################

        numVars = len(outs)
        for i in range(0, numVars):
            dendro.generate_separate([outs[i]], [vnames[i]], "[pp]")


def main(staged_type, gauge, eta_damp, prefix, enable_ssl, enable_cahd):
    if enable_ssl:
        print("// CODEGEN: SSL was enabled, adding term to gauge condition!")

    if enable_cahd:
        print("// CODEGEN: CAHD was enabled, adding damping term to chi!")

    if staged_type == "staged":
        print("//Codgen: generating staged version ")
        if gauge == "rochester":
            print("//Codgen: using rochester gauge")
            if eta_damp == "func":
                print("//Codgen: using eta func damping")
                bssn_rochester_puncture_gauge(eta_func, True, prefix, enable_ssl)
            else:
                print("//Codgen: using eta const damping")
                bssn_rochester_puncture_gauge(eta, True, prefix, enable_ssl)

        else:
            print("//Codgen: using standard gauge")
            if eta_damp == "func":
                print("//Codgen: using eta func damping")
                bssn_puncture_gauge(eta_func, True, prefix, enable_ssl, enable_cahd)
            else:
                print("//Codgen: using eta const damping")
                bssn_puncture_gauge(eta, True, prefix, enable_ssl, enable_cahd)

    else:
        print("//Codgen: generating unstage version ")
        if gauge == "rochester":
            print("//Codgen: using rochester gauge")
            if eta_damp == "func":
                print("//Codgen: using eta func damping")
                bssn_rochester_puncture_gauge(eta_func, False, prefix, enable_ssl)
            else:
                print("//Codgen: using eta const damping")
                bssn_rochester_puncture_gauge(eta, False, prefix, enable_ssl)

        else:
            print("//Codgen: using standard gauge")
            if eta_damp == "func":
                print("//Codgen: using eta func damping")
                bssn_puncture_gauge(eta_func, False, prefix, enable_ssl, enable_cahd)
            else:
                print("//Codgen: using eta const damping")
                bssn_puncture_gauge(eta, False, prefix, enable_ssl, enable_cahd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="BSSN Code Generation",
        description="Generate the Code for the BSSN RHS equations.",
    )

    parser.add_argument(
        "-t",
        "--staged_type",
        choices=["staged", "unstaged"],
        default="unstaged",
        help="If we should use staged or unstaged code",
    )
    parser.add_argument(
        "-g",
        "--gauge",
        choices=["standard", "rochester"],
        default="standard",
        help="The gauge type",
    )
    parser.add_argument(
        "-e",
        "--eta_damp",
        choices=["const", "func"],
        default="const",
        help="The eta damping type, a function or a constant",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        help="The file prefix for staged version",
        default="output_",
    )
    parser.add_argument(
        "-s",
        "--enable_ssl",
        action="store_true",
        help="Whether or not to generate with the SSL code",
    )
    parser.add_argument(
        "-c",
        "--enable_cahd",
        action="store_true",
        help="Whether or not to enable CAHD",
    )

    args = parser.parse_args()

    main(**vars(args))

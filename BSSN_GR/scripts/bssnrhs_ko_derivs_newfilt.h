if (bssn::BSSN_DERIVS->get_filter_family() ==
    dendroderivs::FilterFamily::FF_KO) {
    // for KO diss, input/output should be the same!
    bssn::BSSN_DERIVS->filter(alpha, a_rhs, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);

    bssn::BSSN_DERIVS->filter(chi, chi_rhs, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);

    bssn::BSSN_DERIVS->filter(K, K_rhs, grad_0_gt0, grad_1_gt0, grad_2_gt0, hx,
                              hy, hz, sigma, sz, bflag);

    bssn::BSSN_DERIVS->filter(gt0, gt_rhs00, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(gt1, gt_rhs01, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(gt2, gt_rhs02, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(gt3, gt_rhs11, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(gt4, gt_rhs12, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(gt5, gt_rhs22, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);

    bssn::BSSN_DERIVS->filter(beta0, b_rhs0, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(beta1, b_rhs1, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(beta2, b_rhs2, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);

    bssn::BSSN_DERIVS->filter(At0, At_rhs00, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(At1, At_rhs01, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(At2, At_rhs02, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(At3, At_rhs11, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(At4, At_rhs12, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(At5, At_rhs22, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);

    bssn::BSSN_DERIVS->filter(Gt0, Gt_rhs0, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(Gt1, Gt_rhs1, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(Gt2, Gt_rhs2, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);

    bssn::BSSN_DERIVS->filter(B0, B_rhs0, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(B1, B_rhs1, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(B2, B_rhs2, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);
}

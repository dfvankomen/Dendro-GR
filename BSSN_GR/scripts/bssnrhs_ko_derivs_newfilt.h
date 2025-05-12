if (bssn::BSSN_DERIVS->get_filter_family() ==
    dendroderivs::FilterFamily::FF_KO) {
    // for KO diss, input/output should be the same!
    bssn::BSSN_DERIVS->filter(a_rhs, a_rhs, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);

    bssn::BSSN_DERIVS->filter(chi_rhs, chi_rhs, grad_0_gt0, grad_1_gt0,
                              grad_2_gt0, hx, hy, hz, sigma, sz, bflag);

    bssn::BSSN_DERIVS->filter(K_rhs, K_rhs, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                              hx, hy, hz, sigma, sz, bflag);

    bssn::BSSN_DERIVS->filter(gt_rhs00, gt_rhs00, grad_0_gt0, grad_1_gt0,
                              grad_2_gt0, hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(gt_rhs01, gt_rhs01, grad_0_gt0, grad_1_gt0,
                              grad_2_gt0, hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(gt_rhs02, gt_rhs02, grad_0_gt0, grad_1_gt0,
                              grad_2_gt0, hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(gt_rhs11, gt_rhs11, grad_0_gt0, grad_1_gt0,
                              grad_2_gt0, hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(gt_rhs12, gt_rhs12, grad_0_gt0, grad_1_gt0,
                              grad_2_gt0, hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(gt_rhs22, gt_rhs22, grad_0_gt0, grad_1_gt0,
                              grad_2_gt0, hx, hy, hz, sigma, sz, bflag);

    bssn::BSSN_DERIVS->filter(b_rhs0, b_rhs0, grad_0_gt0, grad_1_gt0,
                              grad_2_gt0, hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(b_rhs1, b_rhs1, grad_0_gt0, grad_1_gt0,
                              grad_2_gt0, hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(b_rhs2, b_rhs2, grad_0_gt0, grad_1_gt0,
                              grad_2_gt0, hx, hy, hz, sigma, sz, bflag);

    bssn::BSSN_DERIVS->filter(At_rhs00, At_rhs00, grad_0_gt0, grad_1_gt0,
                              grad_2_gt0, hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(At_rhs01, At_rhs01, grad_0_gt0, grad_1_gt0,
                              grad_2_gt0, hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(At_rhs02, At_rhs02, grad_0_gt0, grad_1_gt0,
                              grad_2_gt0, hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(At_rhs11, At_rhs11, grad_0_gt0, grad_1_gt0,
                              grad_2_gt0, hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(At_rhs12, At_rhs12, grad_0_gt0, grad_1_gt0,
                              grad_2_gt0, hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(At_rhs22, At_rhs22, grad_0_gt0, grad_1_gt0,
                              grad_2_gt0, hx, hy, hz, sigma, sz, bflag);

    bssn::BSSN_DERIVS->filter(Gt_rhs0, Gt_rhs0, grad_0_gt0, grad_1_gt0,
                              grad_2_gt0, hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(Gt_rhs1, Gt_rhs1, grad_0_gt0, grad_1_gt0,
                              grad_2_gt0, hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(Gt_rhs2, Gt_rhs2, grad_0_gt0, grad_1_gt0,
                              grad_2_gt0, hx, hy, hz, sigma, sz, bflag);

    bssn::BSSN_DERIVS->filter(B_rhs0, B_rhs0, grad_0_gt0, grad_1_gt0,
                              grad_2_gt0, hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(B_rhs1, B_rhs1, grad_0_gt0, grad_1_gt0,
                              grad_2_gt0, hx, hy, hz, sigma, sz, bflag);
    bssn::BSSN_DERIVS->filter(B_rhs2, B_rhs2, grad_0_gt0, grad_1_gt0,
                              grad_2_gt0, hx, hy, hz, sigma, sz, bflag);
}

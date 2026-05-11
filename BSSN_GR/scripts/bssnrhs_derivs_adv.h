#ifdef BSSN_USE_ADVECTIVE_DERIVS
adv_deriv_x(agrad_0_gt0, gt0, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_gt0, gt0, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_gt0, gt0, hz, sz, beta2, bflag);
adv_deriv_x(agrad_0_gt1, gt1, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_gt1, gt1, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_gt1, gt1, hz, sz, beta2, bflag);
adv_deriv_x(agrad_0_gt2, gt2, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_gt2, gt2, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_gt2, gt2, hz, sz, beta2, bflag);
adv_deriv_x(agrad_0_gt3, gt3, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_gt3, gt3, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_gt3, gt3, hz, sz, beta2, bflag);
adv_deriv_x(agrad_0_gt4, gt4, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_gt4, gt4, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_gt4, gt4, hz, sz, beta2, bflag);
adv_deriv_x(agrad_0_gt5, gt5, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_gt5, gt5, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_gt5, gt5, hz, sz, beta2, bflag);
adv_deriv_x(agrad_0_At0, At0, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_At0, At0, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_At0, At0, hz, sz, beta2, bflag);
adv_deriv_x(agrad_0_At1, At1, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_At1, At1, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_At1, At1, hz, sz, beta2, bflag);
adv_deriv_x(agrad_0_At2, At2, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_At2, At2, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_At2, At2, hz, sz, beta2, bflag);
adv_deriv_x(agrad_0_At3, At3, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_At3, At3, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_At3, At3, hz, sz, beta2, bflag);
adv_deriv_x(agrad_0_At4, At4, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_At4, At4, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_At4, At4, hz, sz, beta2, bflag);
adv_deriv_x(agrad_0_At5, At5, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_At5, At5, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_At5, At5, hz, sz, beta2, bflag);
adv_deriv_x(agrad_0_alpha, alpha, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_alpha, alpha, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_alpha, alpha, hz, sz, beta2, bflag);
adv_deriv_x(agrad_0_beta0, beta0, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_beta0, beta0, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_beta0, beta0, hz, sz, beta2, bflag);
adv_deriv_x(agrad_0_beta1, beta1, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_beta1, beta1, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_beta1, beta1, hz, sz, beta2, bflag);
adv_deriv_x(agrad_0_beta2, beta2, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_beta2, beta2, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_beta2, beta2, hz, sz, beta2, bflag);
adv_deriv_x(agrad_0_chi, chi, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_chi, chi, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_chi, chi, hz, sz, beta2, bflag);
adv_deriv_x(agrad_0_Gt0, Gt0, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_Gt0, Gt0, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_Gt0, Gt0, hz, sz, beta2, bflag);
adv_deriv_x(agrad_0_Gt1, Gt1, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_Gt1, Gt1, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_Gt1, Gt1, hz, sz, beta2, bflag);
adv_deriv_x(agrad_0_Gt2, Gt2, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_Gt2, Gt2, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_Gt2, Gt2, hz, sz, beta2, bflag);
double* agrad_0_Gh0 = agrad_0_Gt0;
double* agrad_1_Gh0 = agrad_1_Gt0;
double* agrad_2_Gh0 = agrad_2_Gt0;
double* agrad_0_Gh1 = agrad_0_Gt1;
double* agrad_1_Gh1 = agrad_1_Gt1;
double* agrad_2_Gh1 = agrad_2_Gt1;
double* agrad_0_Gh2 = agrad_0_Gt2;
double* agrad_1_Gh2 = agrad_1_Gt2;
double* agrad_2_Gh2 = agrad_2_Gt2;
double* agrad_0_Gammahat0 = agrad_0_Gt0;
double* agrad_1_Gammahat0 = agrad_1_Gt0;
double* agrad_2_Gammahat0 = agrad_2_Gt0;
double* agrad_0_Gammahat1 = agrad_0_Gt1;
double* agrad_1_Gammahat1 = agrad_1_Gt1;
double* agrad_2_Gammahat1 = agrad_2_Gt1;
double* agrad_0_Gammahat2 = agrad_0_Gt2;
double* agrad_1_Gammahat2 = agrad_1_Gt2;
double* agrad_2_Gammahat2 = agrad_2_Gt2;
adv_deriv_x(agrad_0_K, K, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_K, K, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_K, K, hz, sz, beta2, bflag);
adv_deriv_x(agrad_0_B0, B0, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_B0, B0, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_B0, B0, hz, sz, beta2, bflag);
adv_deriv_x(agrad_0_B1, B1, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_B1, B1, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_B1, B1, hz, sz, beta2, bflag);
adv_deriv_x(agrad_0_B2, B2, hx, sz, beta0, bflag);
adv_deriv_y(agrad_1_B2, B2, hy, sz, beta1, bflag);
adv_deriv_z(agrad_2_B2, B2, hz, sz, beta2, bflag);

#else

// Fallback: use centered derivatives when advective derivatives are disabled.
// This keeps generated CCZ4 code compiling because it always references agrad_*.

double* agrad_0_gt0 = grad_0_gt0;
double* agrad_1_gt0 = grad_1_gt0;
double* agrad_2_gt0 = grad_2_gt0;
double* agrad_0_gt1 = grad_0_gt1;
double* agrad_1_gt1 = grad_1_gt1;
double* agrad_2_gt1 = grad_2_gt1;
double* agrad_0_gt2 = grad_0_gt2;
double* agrad_1_gt2 = grad_1_gt2;
double* agrad_2_gt2 = grad_2_gt2;
double* agrad_0_gt3 = grad_0_gt3;
double* agrad_1_gt3 = grad_1_gt3;
double* agrad_2_gt3 = grad_2_gt3;
double* agrad_0_gt4 = grad_0_gt4;
double* agrad_1_gt4 = grad_1_gt4;
double* agrad_2_gt4 = grad_2_gt4;
double* agrad_0_gt5 = grad_0_gt5;
double* agrad_1_gt5 = grad_1_gt5;
double* agrad_2_gt5 = grad_2_gt5;

double* agrad_0_At0 = grad_0_At0;
double* agrad_1_At0 = grad_1_At0;
double* agrad_2_At0 = grad_2_At0;
double* agrad_0_At1 = grad_0_At1;
double* agrad_1_At1 = grad_1_At1;
double* agrad_2_At1 = grad_2_At1;
double* agrad_0_At2 = grad_0_At2;
double* agrad_1_At2 = grad_1_At2;
double* agrad_2_At2 = grad_2_At2;
double* agrad_0_At3 = grad_0_At3;
double* agrad_1_At3 = grad_1_At3;
double* agrad_2_At3 = grad_2_At3;
double* agrad_0_At4 = grad_0_At4;
double* agrad_1_At4 = grad_1_At4;
double* agrad_2_At4 = grad_2_At4;
double* agrad_0_At5 = grad_0_At5;
double* agrad_1_At5 = grad_1_At5;
double* agrad_2_At5 = grad_2_At5;

double* agrad_0_alpha = grad_0_alpha;
double* agrad_1_alpha = grad_1_alpha;
double* agrad_2_alpha = grad_2_alpha;

double* agrad_0_beta0 = grad_0_beta0;
double* agrad_1_beta0 = grad_1_beta0;
double* agrad_2_beta0 = grad_2_beta0;
double* agrad_0_beta1 = grad_0_beta1;
double* agrad_1_beta1 = grad_1_beta1;
double* agrad_2_beta1 = grad_2_beta1;
double* agrad_0_beta2 = grad_0_beta2;
double* agrad_1_beta2 = grad_1_beta2;
double* agrad_2_beta2 = grad_2_beta2;

double* agrad_0_chi = grad_0_chi;
double* agrad_1_chi = grad_1_chi;
double* agrad_2_chi = grad_2_chi;

double* agrad_0_K = grad_0_K;
double* agrad_1_K = grad_1_K;
double* agrad_2_K = grad_2_K;

double* agrad_0_B0 = grad_0_B0;
double* agrad_1_B0 = grad_1_B0;
double* agrad_2_B0 = grad_2_B0;
double* agrad_0_B1 = grad_0_B1;
double* agrad_1_B1 = grad_1_B1;
double* agrad_2_B1 = grad_2_B1;
double* agrad_0_B2 = grad_0_B2;
double* agrad_1_B2 = grad_1_B2;
double* agrad_2_B2 = grad_2_B2;

// Your code seems to use Gt as the actual allocated array name.
// Keep all aliases consistent.
double* agrad_0_Gt0 = grad_0_Gt0;
double* agrad_1_Gt0 = grad_1_Gt0;
double* agrad_2_Gt0 = grad_2_Gt0;
double* agrad_0_Gt1 = grad_0_Gt1;
double* agrad_1_Gt1 = grad_1_Gt1;
double* agrad_2_Gt1 = grad_2_Gt1;
double* agrad_0_Gt2 = grad_0_Gt2;
double* agrad_1_Gt2 = grad_1_Gt2;
double* agrad_2_Gt2 = grad_2_Gt2;

double* agrad_0_Gh0 = agrad_0_Gt0;
double* agrad_1_Gh0 = agrad_1_Gt0;
double* agrad_2_Gh0 = agrad_2_Gt0;
double* agrad_0_Gh1 = agrad_0_Gt1;
double* agrad_1_Gh1 = agrad_1_Gt1;
double* agrad_2_Gh1 = agrad_2_Gt1;
double* agrad_0_Gh2 = agrad_0_Gt2;
double* agrad_1_Gh2 = agrad_1_Gt2;
double* agrad_2_Gh2 = agrad_2_Gt2;

double* agrad_0_Gammahat0 = agrad_0_Gt0;
double* agrad_1_Gammahat0 = agrad_1_Gt0;
double* agrad_2_Gammahat0 = agrad_2_Gt0;
double* agrad_0_Gammahat1 = agrad_0_Gt1;
double* agrad_1_Gammahat1 = agrad_1_Gt1;
double* agrad_2_Gammahat1 = agrad_2_Gt1;
double* agrad_0_Gammahat2 = agrad_0_Gt2;
double* agrad_1_Gammahat2 = agrad_1_Gt2;
double* agrad_2_Gammahat2 = agrad_2_Gt2;

#endif
const DEVICE_REAL grad_0_alpha =
    *(deriv_evars->grad_0_alpha + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_alpha =
    *(deriv_evars->grad_1_alpha + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_alpha =
    *(deriv_evars->grad_2_alpha + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_chi =
    *(deriv_evars->grad_0_chi + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_chi =
    *(deriv_evars->grad_1_chi + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_chi =
    *(deriv_evars->grad_2_chi + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_K = *(deriv_evars->grad_0_K + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_K = *(deriv_evars->grad_1_K + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_K = *(deriv_evars->grad_2_K + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_Gt0 =
    *(deriv_evars->grad_0_Gt0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_Gt0 =
    *(deriv_evars->grad_1_Gt0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_Gt0 =
    *(deriv_evars->grad_2_Gt0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_Gt1 =
    *(deriv_evars->grad_0_Gt1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_Gt1 =
    *(deriv_evars->grad_1_Gt1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_Gt1 =
    *(deriv_evars->grad_2_Gt1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_Gt2 =
    *(deriv_evars->grad_0_Gt2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_Gt2 =
    *(deriv_evars->grad_1_Gt2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_Gt2 =
    *(deriv_evars->grad_2_Gt2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_beta0 =
    *(deriv_evars->grad_0_beta0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_beta0 =
    *(deriv_evars->grad_1_beta0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_beta0 =
    *(deriv_evars->grad_2_beta0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_beta1 =
    *(deriv_evars->grad_0_beta1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_beta1 =
    *(deriv_evars->grad_1_beta1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_beta1 =
    *(deriv_evars->grad_2_beta1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_beta2 =
    *(deriv_evars->grad_0_beta2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_beta2 =
    *(deriv_evars->grad_1_beta2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_beta2 =
    *(deriv_evars->grad_2_beta2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_B0 = *(deriv_evars->grad_0_B0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_B0 = *(deriv_evars->grad_1_B0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_B0 = *(deriv_evars->grad_2_B0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_B1 = *(deriv_evars->grad_0_B1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_B1 = *(deriv_evars->grad_1_B1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_B1 = *(deriv_evars->grad_2_B1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_B2 = *(deriv_evars->grad_0_B2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_B2 = *(deriv_evars->grad_1_B2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_B2 = *(deriv_evars->grad_2_B2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_gt0 =
    *(deriv_evars->grad_0_gt0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_gt0 =
    *(deriv_evars->grad_1_gt0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_gt0 =
    *(deriv_evars->grad_2_gt0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_gt1 =
    *(deriv_evars->grad_0_gt1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_gt1 =
    *(deriv_evars->grad_1_gt1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_gt1 =
    *(deriv_evars->grad_2_gt1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_gt2 =
    *(deriv_evars->grad_0_gt2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_gt2 =
    *(deriv_evars->grad_1_gt2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_gt2 =
    *(deriv_evars->grad_2_gt2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_gt3 =
    *(deriv_evars->grad_0_gt3 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_gt3 =
    *(deriv_evars->grad_1_gt3 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_gt3 =
    *(deriv_evars->grad_2_gt3 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_gt4 =
    *(deriv_evars->grad_0_gt4 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_gt4 =
    *(deriv_evars->grad_1_gt4 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_gt4 =
    *(deriv_evars->grad_2_gt4 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_gt5 =
    *(deriv_evars->grad_0_gt5 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_gt5 =
    *(deriv_evars->grad_1_gt5 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_gt5 =
    *(deriv_evars->grad_2_gt5 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_At0 =
    *(deriv_evars->grad_0_At0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_At0 =
    *(deriv_evars->grad_1_At0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_At0 =
    *(deriv_evars->grad_2_At0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_At1 =
    *(deriv_evars->grad_0_At1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_At1 =
    *(deriv_evars->grad_1_At1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_At1 =
    *(deriv_evars->grad_2_At1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_At2 =
    *(deriv_evars->grad_0_At2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_At2 =
    *(deriv_evars->grad_1_At2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_At2 =
    *(deriv_evars->grad_2_At2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_At3 =
    *(deriv_evars->grad_0_At3 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_At3 =
    *(deriv_evars->grad_1_At3 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_At3 =
    *(deriv_evars->grad_2_At3 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_At4 =
    *(deriv_evars->grad_0_At4 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_At4 =
    *(deriv_evars->grad_1_At4 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_At4 =
    *(deriv_evars->grad_2_At4 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_0_At5 =
    *(deriv_evars->grad_0_At5 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_1_At5 =
    *(deriv_evars->grad_1_At5 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad_2_At5 =
    *(deriv_evars->grad_2_At5 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_0_gt0 =
    *(deriv_evars->grad2_0_0_gt0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_1_gt0 =
    *(deriv_evars->grad2_0_1_gt0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_2_gt0 =
    *(deriv_evars->grad2_0_2_gt0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_1_1_gt0 =
    *(deriv_evars->grad2_1_1_gt0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_1_2_gt0 =
    *(deriv_evars->grad2_1_2_gt0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_2_2_gt0 =
    *(deriv_evars->grad2_2_2_gt0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_0_gt1 =
    *(deriv_evars->grad2_0_0_gt1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_1_gt1 =
    *(deriv_evars->grad2_0_1_gt1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_2_gt1 =
    *(deriv_evars->grad2_0_2_gt1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_1_1_gt1 =
    *(deriv_evars->grad2_1_1_gt1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_1_2_gt1 =
    *(deriv_evars->grad2_1_2_gt1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_2_2_gt1 =
    *(deriv_evars->grad2_2_2_gt1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_0_gt2 =
    *(deriv_evars->grad2_0_0_gt2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_1_gt2 =
    *(deriv_evars->grad2_0_1_gt2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_2_gt2 =
    *(deriv_evars->grad2_0_2_gt2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_1_1_gt2 =
    *(deriv_evars->grad2_1_1_gt2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_1_2_gt2 =
    *(deriv_evars->grad2_1_2_gt2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_2_2_gt2 =
    *(deriv_evars->grad2_2_2_gt2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_0_gt3 =
    *(deriv_evars->grad2_0_0_gt3 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_1_gt3 =
    *(deriv_evars->grad2_0_1_gt3 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_2_gt3 =
    *(deriv_evars->grad2_0_2_gt3 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_1_1_gt3 =
    *(deriv_evars->grad2_1_1_gt3 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_1_2_gt3 =
    *(deriv_evars->grad2_1_2_gt3 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_2_2_gt3 =
    *(deriv_evars->grad2_2_2_gt3 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_0_gt4 =
    *(deriv_evars->grad2_0_0_gt4 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_1_gt4 =
    *(deriv_evars->grad2_0_1_gt4 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_2_gt4 =
    *(deriv_evars->grad2_0_2_gt4 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_1_1_gt4 =
    *(deriv_evars->grad2_1_1_gt4 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_1_2_gt4 =
    *(deriv_evars->grad2_1_2_gt4 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_2_2_gt4 =
    *(deriv_evars->grad2_2_2_gt4 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_0_gt5 =
    *(deriv_evars->grad2_0_0_gt5 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_1_gt5 =
    *(deriv_evars->grad2_0_1_gt5 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_2_gt5 =
    *(deriv_evars->grad2_0_2_gt5 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_1_1_gt5 =
    *(deriv_evars->grad2_1_1_gt5 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_1_2_gt5 =
    *(deriv_evars->grad2_1_2_gt5 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_2_2_gt5 =
    *(deriv_evars->grad2_2_2_gt5 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_0_chi =
    *(deriv_evars->grad2_0_0_chi + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_1_chi =
    *(deriv_evars->grad2_0_1_chi + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_2_chi =
    *(deriv_evars->grad2_0_2_chi + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_1_1_chi =
    *(deriv_evars->grad2_1_1_chi + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_1_2_chi =
    *(deriv_evars->grad2_1_2_chi + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_2_2_chi =
    *(deriv_evars->grad2_2_2_chi + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_0_alpha =
    *(deriv_evars->grad2_0_0_alpha + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_1_alpha =
    *(deriv_evars->grad2_0_1_alpha + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_2_alpha =
    *(deriv_evars->grad2_0_2_alpha + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_1_1_alpha =
    *(deriv_evars->grad2_1_1_alpha + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_1_2_alpha =
    *(deriv_evars->grad2_1_2_alpha + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_2_2_alpha =
    *(deriv_evars->grad2_2_2_alpha + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_0_beta0 =
    *(deriv_evars->grad2_0_0_beta0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_1_beta0 =
    *(deriv_evars->grad2_0_1_beta0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_2_beta0 =
    *(deriv_evars->grad2_0_2_beta0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_1_1_beta0 =
    *(deriv_evars->grad2_1_1_beta0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_1_2_beta0 =
    *(deriv_evars->grad2_1_2_beta0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_2_2_beta0 =
    *(deriv_evars->grad2_2_2_beta0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_0_beta1 =
    *(deriv_evars->grad2_0_0_beta1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_1_beta1 =
    *(deriv_evars->grad2_0_1_beta1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_2_beta1 =
    *(deriv_evars->grad2_0_2_beta1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_1_1_beta1 =
    *(deriv_evars->grad2_1_1_beta1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_1_2_beta1 =
    *(deriv_evars->grad2_1_2_beta1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_2_2_beta1 =
    *(deriv_evars->grad2_2_2_beta1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_0_beta2 =
    *(deriv_evars->grad2_0_0_beta2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_1_beta2 =
    *(deriv_evars->grad2_0_1_beta2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_0_2_beta2 =
    *(deriv_evars->grad2_0_2_beta2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_1_1_beta2 =
    *(deriv_evars->grad2_1_1_beta2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_1_2_beta2 =
    *(deriv_evars->grad2_1_2_beta2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL grad2_2_2_beta2 =
    *(deriv_evars->grad2_2_2_beta2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_alpha =
    *(deriv_evars->kograd_0_alpha + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_alpha =
    *(deriv_evars->kograd_1_alpha + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_alpha =
    *(deriv_evars->kograd_2_alpha + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_chi =
    *(deriv_evars->kograd_0_chi + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_chi =
    *(deriv_evars->kograd_1_chi + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_chi =
    *(deriv_evars->kograd_2_chi + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_K =
    *(deriv_evars->kograd_0_K + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_K =
    *(deriv_evars->kograd_1_K + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_K =
    *(deriv_evars->kograd_2_K + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_Gt0 =
    *(deriv_evars->kograd_0_Gt0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_Gt0 =
    *(deriv_evars->kograd_1_Gt0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_Gt0 =
    *(deriv_evars->kograd_2_Gt0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_Gt1 =
    *(deriv_evars->kograd_0_Gt1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_Gt1 =
    *(deriv_evars->kograd_1_Gt1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_Gt1 =
    *(deriv_evars->kograd_2_Gt1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_Gt2 =
    *(deriv_evars->kograd_0_Gt2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_Gt2 =
    *(deriv_evars->kograd_1_Gt2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_Gt2 =
    *(deriv_evars->kograd_2_Gt2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_beta0 =
    *(deriv_evars->kograd_0_beta0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_beta0 =
    *(deriv_evars->kograd_1_beta0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_beta0 =
    *(deriv_evars->kograd_2_beta0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_beta1 =
    *(deriv_evars->kograd_0_beta1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_beta1 =
    *(deriv_evars->kograd_1_beta1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_beta1 =
    *(deriv_evars->kograd_2_beta1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_beta2 =
    *(deriv_evars->kograd_0_beta2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_beta2 =
    *(deriv_evars->kograd_1_beta2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_beta2 =
    *(deriv_evars->kograd_2_beta2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_B0 =
    *(deriv_evars->kograd_0_B0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_B0 =
    *(deriv_evars->kograd_1_B0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_B0 =
    *(deriv_evars->kograd_2_B0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_B1 =
    *(deriv_evars->kograd_0_B1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_B1 =
    *(deriv_evars->kograd_1_B1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_B1 =
    *(deriv_evars->kograd_2_B1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_B2 =
    *(deriv_evars->kograd_0_B2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_B2 =
    *(deriv_evars->kograd_1_B2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_B2 =
    *(deriv_evars->kograd_2_B2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_gt0 =
    *(deriv_evars->kograd_0_gt0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_gt0 =
    *(deriv_evars->kograd_1_gt0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_gt0 =
    *(deriv_evars->kograd_2_gt0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_gt1 =
    *(deriv_evars->kograd_0_gt1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_gt1 =
    *(deriv_evars->kograd_1_gt1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_gt1 =
    *(deriv_evars->kograd_2_gt1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_gt2 =
    *(deriv_evars->kograd_0_gt2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_gt2 =
    *(deriv_evars->kograd_1_gt2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_gt2 =
    *(deriv_evars->kograd_2_gt2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_gt3 =
    *(deriv_evars->kograd_0_gt3 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_gt3 =
    *(deriv_evars->kograd_1_gt3 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_gt3 =
    *(deriv_evars->kograd_2_gt3 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_gt4 =
    *(deriv_evars->kograd_0_gt4 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_gt4 =
    *(deriv_evars->kograd_1_gt4 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_gt4 =
    *(deriv_evars->kograd_2_gt4 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_gt5 =
    *(deriv_evars->kograd_0_gt5 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_gt5 =
    *(deriv_evars->kograd_1_gt5 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_gt5 =
    *(deriv_evars->kograd_2_gt5 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_At0 =
    *(deriv_evars->kograd_0_At0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_At0 =
    *(deriv_evars->kograd_1_At0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_At0 =
    *(deriv_evars->kograd_2_At0 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_At1 =
    *(deriv_evars->kograd_0_At1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_At1 =
    *(deriv_evars->kograd_1_At1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_At1 =
    *(deriv_evars->kograd_2_At1 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_At2 =
    *(deriv_evars->kograd_0_At2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_At2 =
    *(deriv_evars->kograd_1_At2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_At2 =
    *(deriv_evars->kograd_2_At2 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_At3 =
    *(deriv_evars->kograd_0_At3 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_At3 =
    *(deriv_evars->kograd_1_At3 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_At3 =
    *(deriv_evars->kograd_2_At3 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_At4 =
    *(deriv_evars->kograd_0_At4 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_At4 =
    *(deriv_evars->kograd_1_At4 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_At4 =
    *(deriv_evars->kograd_2_At4 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_0_At5 =
    *(deriv_evars->kograd_0_At5 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_1_At5 =
    *(deriv_evars->kograd_1_At5 + BLK_ID * BLK_SZ + pp);
const DEVICE_REAL kograd_2_At5 =
    *(deriv_evars->kograd_2_At5 + BLK_ID * BLK_SZ + pp);

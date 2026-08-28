// hand-generated from bssnrhs_derivs.h: one planned grad_set per variable
// (dendrolib DendroDerivatives::grad_set); the engine picks intermediate vs
// terminal shapes. Interior blocks only; puncture blocks keep bssnrhs_derivs.h.
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_alpha;
    o.y = grad_1_alpha;
    o.z = grad_2_alpha;
    o.xx = grad2_0_0_alpha;
    o.yy = grad2_1_1_alpha;
    o.zz = grad2_2_2_alpha;
    o.xy = grad2_0_1_alpha;
    o.xz = grad2_0_2_alpha;
    o.yz = grad2_1_2_alpha;
    bssn::active_derivs()->grad_set(o, alpha, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z | dendroderivs::DendroDerivatives::DM_XX | dendroderivs::DendroDerivatives::DM_YY | dendroderivs::DendroDerivatives::DM_ZZ | dendroderivs::DendroDerivatives::DM_XY | dendroderivs::DendroDerivatives::DM_XZ | dendroderivs::DendroDerivatives::DM_YZ, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_beta0;
    o.y = grad_1_beta0;
    o.z = grad_2_beta0;
    o.xx = grad2_0_0_beta0;
    o.yy = grad2_1_1_beta0;
    o.zz = grad2_2_2_beta0;
    o.xy = grad2_0_1_beta0;
    o.xz = grad2_0_2_beta0;
    o.yz = grad2_1_2_beta0;
    bssn::active_derivs()->grad_set(o, beta0, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z | dendroderivs::DendroDerivatives::DM_XX | dendroderivs::DendroDerivatives::DM_YY | dendroderivs::DendroDerivatives::DM_ZZ | dendroderivs::DendroDerivatives::DM_XY | dendroderivs::DendroDerivatives::DM_XZ | dendroderivs::DendroDerivatives::DM_YZ, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_beta1;
    o.y = grad_1_beta1;
    o.z = grad_2_beta1;
    o.xx = grad2_0_0_beta1;
    o.yy = grad2_1_1_beta1;
    o.zz = grad2_2_2_beta1;
    o.xy = grad2_0_1_beta1;
    o.xz = grad2_0_2_beta1;
    o.yz = grad2_1_2_beta1;
    bssn::active_derivs()->grad_set(o, beta1, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z | dendroderivs::DendroDerivatives::DM_XX | dendroderivs::DendroDerivatives::DM_YY | dendroderivs::DendroDerivatives::DM_ZZ | dendroderivs::DendroDerivatives::DM_XY | dendroderivs::DendroDerivatives::DM_XZ | dendroderivs::DendroDerivatives::DM_YZ, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_beta2;
    o.y = grad_1_beta2;
    o.z = grad_2_beta2;
    o.xx = grad2_0_0_beta2;
    o.yy = grad2_1_1_beta2;
    o.zz = grad2_2_2_beta2;
    o.xy = grad2_0_1_beta2;
    o.xz = grad2_0_2_beta2;
    o.yz = grad2_1_2_beta2;
    bssn::active_derivs()->grad_set(o, beta2, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z | dendroderivs::DendroDerivatives::DM_XX | dendroderivs::DendroDerivatives::DM_YY | dendroderivs::DendroDerivatives::DM_ZZ | dendroderivs::DendroDerivatives::DM_XY | dendroderivs::DendroDerivatives::DM_XZ | dendroderivs::DendroDerivatives::DM_YZ, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_B0;
    o.y = grad_1_B0;
    o.z = grad_2_B0;
    bssn::active_derivs()->grad_set(o, B0, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_B1;
    o.y = grad_1_B1;
    o.z = grad_2_B1;
    bssn::active_derivs()->grad_set(o, B1, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_B2;
    o.y = grad_1_B2;
    o.z = grad_2_B2;
    bssn::active_derivs()->grad_set(o, B2, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_chi;
    o.y = grad_1_chi;
    o.z = grad_2_chi;
    o.xx = grad2_0_0_chi;
    o.yy = grad2_1_1_chi;
    o.zz = grad2_2_2_chi;
    o.xy = grad2_0_1_chi;
    o.xz = grad2_0_2_chi;
    o.yz = grad2_1_2_chi;
    bssn::active_derivs()->grad_set(o, chi, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z | dendroderivs::DendroDerivatives::DM_XX | dendroderivs::DendroDerivatives::DM_YY | dendroderivs::DendroDerivatives::DM_ZZ | dendroderivs::DendroDerivatives::DM_XY | dendroderivs::DendroDerivatives::DM_XZ | dendroderivs::DendroDerivatives::DM_YZ, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_Gt0;
    o.y = grad_1_Gt0;
    o.z = grad_2_Gt0;
    bssn::active_derivs()->grad_set(o, Gt0, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_Gt1;
    o.y = grad_1_Gt1;
    o.z = grad_2_Gt1;
    bssn::active_derivs()->grad_set(o, Gt1, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_Gt2;
    o.y = grad_1_Gt2;
    o.z = grad_2_Gt2;
    bssn::active_derivs()->grad_set(o, Gt2, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_K;
    o.y = grad_1_K;
    o.z = grad_2_K;
    bssn::active_derivs()->grad_set(o, K, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_gt0;
    o.y = grad_1_gt0;
    o.z = grad_2_gt0;
    o.xx = grad2_0_0_gt0;
    o.yy = grad2_1_1_gt0;
    o.zz = grad2_2_2_gt0;
    o.xy = grad2_0_1_gt0;
    o.xz = grad2_0_2_gt0;
    o.yz = grad2_1_2_gt0;
    bssn::active_derivs()->grad_set(o, gt0, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z | dendroderivs::DendroDerivatives::DM_XX | dendroderivs::DendroDerivatives::DM_YY | dendroderivs::DendroDerivatives::DM_ZZ | dendroderivs::DendroDerivatives::DM_XY | dendroderivs::DendroDerivatives::DM_XZ | dendroderivs::DendroDerivatives::DM_YZ, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_gt1;
    o.y = grad_1_gt1;
    o.z = grad_2_gt1;
    o.xx = grad2_0_0_gt1;
    o.yy = grad2_1_1_gt1;
    o.zz = grad2_2_2_gt1;
    o.xy = grad2_0_1_gt1;
    o.xz = grad2_0_2_gt1;
    o.yz = grad2_1_2_gt1;
    bssn::active_derivs()->grad_set(o, gt1, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z | dendroderivs::DendroDerivatives::DM_XX | dendroderivs::DendroDerivatives::DM_YY | dendroderivs::DendroDerivatives::DM_ZZ | dendroderivs::DendroDerivatives::DM_XY | dendroderivs::DendroDerivatives::DM_XZ | dendroderivs::DendroDerivatives::DM_YZ, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_gt2;
    o.y = grad_1_gt2;
    o.z = grad_2_gt2;
    o.xx = grad2_0_0_gt2;
    o.yy = grad2_1_1_gt2;
    o.zz = grad2_2_2_gt2;
    o.xy = grad2_0_1_gt2;
    o.xz = grad2_0_2_gt2;
    o.yz = grad2_1_2_gt2;
    bssn::active_derivs()->grad_set(o, gt2, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z | dendroderivs::DendroDerivatives::DM_XX | dendroderivs::DendroDerivatives::DM_YY | dendroderivs::DendroDerivatives::DM_ZZ | dendroderivs::DendroDerivatives::DM_XY | dendroderivs::DendroDerivatives::DM_XZ | dendroderivs::DendroDerivatives::DM_YZ, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_gt3;
    o.y = grad_1_gt3;
    o.z = grad_2_gt3;
    o.xx = grad2_0_0_gt3;
    o.yy = grad2_1_1_gt3;
    o.zz = grad2_2_2_gt3;
    o.xy = grad2_0_1_gt3;
    o.xz = grad2_0_2_gt3;
    o.yz = grad2_1_2_gt3;
    bssn::active_derivs()->grad_set(o, gt3, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z | dendroderivs::DendroDerivatives::DM_XX | dendroderivs::DendroDerivatives::DM_YY | dendroderivs::DendroDerivatives::DM_ZZ | dendroderivs::DendroDerivatives::DM_XY | dendroderivs::DendroDerivatives::DM_XZ | dendroderivs::DendroDerivatives::DM_YZ, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_gt4;
    o.y = grad_1_gt4;
    o.z = grad_2_gt4;
    o.xx = grad2_0_0_gt4;
    o.yy = grad2_1_1_gt4;
    o.zz = grad2_2_2_gt4;
    o.xy = grad2_0_1_gt4;
    o.xz = grad2_0_2_gt4;
    o.yz = grad2_1_2_gt4;
    bssn::active_derivs()->grad_set(o, gt4, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z | dendroderivs::DendroDerivatives::DM_XX | dendroderivs::DendroDerivatives::DM_YY | dendroderivs::DendroDerivatives::DM_ZZ | dendroderivs::DendroDerivatives::DM_XY | dendroderivs::DendroDerivatives::DM_XZ | dendroderivs::DendroDerivatives::DM_YZ, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_gt5;
    o.y = grad_1_gt5;
    o.z = grad_2_gt5;
    o.xx = grad2_0_0_gt5;
    o.yy = grad2_1_1_gt5;
    o.zz = grad2_2_2_gt5;
    o.xy = grad2_0_1_gt5;
    o.xz = grad2_0_2_gt5;
    o.yz = grad2_1_2_gt5;
    bssn::active_derivs()->grad_set(o, gt5, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z | dendroderivs::DendroDerivatives::DM_XX | dendroderivs::DendroDerivatives::DM_YY | dendroderivs::DendroDerivatives::DM_ZZ | dendroderivs::DendroDerivatives::DM_XY | dendroderivs::DendroDerivatives::DM_XZ | dendroderivs::DendroDerivatives::DM_YZ, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_At0;
    o.y = grad_1_At0;
    o.z = grad_2_At0;
    bssn::active_derivs()->grad_set(o, At0, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_At1;
    o.y = grad_1_At1;
    o.z = grad_2_At1;
    bssn::active_derivs()->grad_set(o, At1, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_At2;
    o.y = grad_1_At2;
    o.z = grad_2_At2;
    bssn::active_derivs()->grad_set(o, At2, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_At3;
    o.y = grad_1_At3;
    o.z = grad_2_At3;
    bssn::active_derivs()->grad_set(o, At3, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_At4;
    o.y = grad_1_At4;
    o.z = grad_2_At4;
    bssn::active_derivs()->grad_set(o, At4, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z, hx, hy, hz, sz, bflag);
}
{
    dendroderivs::DendroDerivatives::DerivSet o;
    o.x = grad_0_At5;
    o.y = grad_1_At5;
    o.z = grad_2_At5;
    bssn::active_derivs()->grad_set(o, At5, dendroderivs::DendroDerivatives::DM_X | dendroderivs::DendroDerivatives::DM_Y | dendroderivs::DendroDerivatives::DM_Z, hx, hy, hz, sz, bflag);
}

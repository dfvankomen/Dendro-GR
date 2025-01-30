//This is an initial data file design to match the one used by Baumgarte: https://arxiv.org/abs/2303.05530
//reading in data from par file...
double A = bssn::TEUK_AMP;
double lambda = bssn::TEUK_WIDTH; //half width really; width is twice this
double r0 = bssn::TEUK_R_0;
//Defining some of the general coordinates and trig functions:
// double xbar   = xx;
// double ybar   = yy;
// double zbar   = zz;
// double rho_sqrd = xbar*xbar + ybar*ybar;
// double rho = sqrt( rho_sqrd );

// double r_sqrd = rho_sqrd + zbar*zbar;
// double rr = sqrt( r_sqrd );

// double epsilon_sq = 1.0e-8;   
// double r_origin = sqrt( r_sqrd + epsilon_sq );

double pi = 3.1415926535;
double E = 2.7182818284590;
double sin_theta;
double cos_theta;


// std::cout << "starting the loop in initial data ... /n";
//double offset = 0.01627605 / 7.0;
double xbar   = xx;
double ybar   = yy;
double zbar   = zz;

double rho_sqrd      = xbar * xbar + ybar * ybar;
double rho           = sqrt(rho_sqrd);

double r_sqrd        = rho_sqrd + zbar * zbar;
double rr            = sqrt(r_sqrd);

double epsilon_sq    = 1.0e-8;  // maybe we can read in grid info later?
double r_origin      = sqrt(r_sqrd + epsilon_sq);

// Away from the axis of symmetry
if (abs(xbar) > 1.0e-7 || abs(ybar) > 1.0e-7) {
    sin_theta  = rho / rr;
    cos_theta  = zbar / rr;
}
// On the axis of symmetry but away from the origin
else if ((abs(xbar) < 1.0e-7 && abs(ybar) < 1.0e-7) && abs(zbar) > 1.0e-7) {
    sin_theta  = rho / rr;
    cos_theta  = zbar / rr;
}
// At the origin
else if ((abs(xbar) < 1.0e-7 && abs(ybar) < 1.0e-7) && abs(zbar) < 1.0e-7) {
    sin_theta  = 0.0;
    cos_theta  = 0.0;
} else {
    std::cout << " problem in initial data /n";
    std::cout << " x = " << xbar;
    std::cout << " y = " << ybar;
    std::cout << " z = " << zbar;
    std::cout << " rr = " << rr;
    std::cout << " rho = " << rho;
}

// // Away from the axis of symmetry
// if ( abs(xbar) > 1.0e-7  ||  abs(ybar) > 1.0e-7 ) {
// 	sin_theta = rho / rr;
// 	cos_theta = zbar / rr;
// } 
// // On the axis of symmetry but away from the origin
// else if ( ( abs(xbar) < 1.0e-7 && abs(ybar) < 1.0e-7 ) && abs(zbar) > 1.0e-7  ) {
// 	sin_theta = rho / rr;
// 	cos_theta = zbar / rr;
// }
// // At the origin
// else if ( ( abs(xbar) < 1.0e-7 && abs(ybar) < 1.0e-7 ) && abs(zbar) < 1.0e-7 ) {
// 	sin_theta = 0.0;
// 	cos_theta = 0.0;
// } else {
// 	std::cout <<" problem in initial data /n";
// 	std::cout << " x = " << xbar ;
// 	std::cout <<" y = "<< ybar ;
// 	std::cout <<" z = "<< zbar;
// 	std::cout << " rr = "<< rr;
// 	std::cout << " rho = "<< rho;
// }
//std::cout << "starting the loop in initial data ... /n";
//Defining the variables we will need later:
double g_xx;
double g_xy;
double g_xz;
double g_yy;
double g_yz;
double g_zz;
double drdx;
double drdy;
double drdz;
double dthdx;
double dthdy;
double dthdz;
double dphdx;
double dphdy;
double dphdz;
double detg;
double rho_axis;
//The Seed function is given by F = 
//These are given in a mathematica script and converted into C code:


double A2;
double B2;
double C2;
//we will here make use of a taylor expanded expression of these functions as to avoid singularities arising from imperfect cancelationsat the origin  We do this out to 
if (rr < 1.0e-2 ) {
A2 = (-64*A*(4*pow(r0,4) - 12*pow(r0,2)*pow(lambda,2) + 3*pow(lambda,4)))/(pow(E,pow(r0,2)/pow(lambda,2))*pow(lambda,4)) + 
   (64*A*pow(rr,2)*(-8*pow(r0,6) + 60*pow(r0,4)*pow(lambda,2) - 90*pow(r0,2)*pow(lambda,4) + 15*pow(lambda,6)))/(5.*pow(E,pow(r0,2)/pow(lambda,2))*pow(lambda,8)) - 
   (32*A*pow(rr,4)*(16*pow(r0,8) - 224*pow(r0,6)*pow(lambda,2) + 840*pow(r0,4)*pow(lambda,4) - 840*pow(r0,2)*pow(lambda,6) + 105*pow(lambda,8)))/
    (35.*pow(E,pow(r0,2)/pow(lambda,2))*pow(lambda,12));
B2 = (-32*A*(4*pow(r0,4) - 12*pow(r0,2)*pow(lambda,2) + 3*pow(lambda,4)))/(pow(E,pow(r0,2)/pow(lambda,2))*pow(lambda,4)) + 
   (32*A*pow(rr,2)*(-8*pow(r0,6) + 60*pow(r0,4)*pow(lambda,2) - 90*pow(r0,2)*pow(lambda,4) + 15*pow(lambda,6)))/(3.*pow(E,pow(r0,2)/pow(lambda,2))*pow(lambda,8)) - 
   (16*A*pow(rr,4)*(16*pow(r0,8) - 224*pow(r0,6)*pow(lambda,2) + 840*pow(r0,4)*pow(lambda,4) - 840*pow(r0,2)*pow(lambda,6) + 105*pow(lambda,8)))/
    (15.*pow(E,pow(r0,2)/pow(lambda,2))*pow(lambda,12));
C2 = (-32*A*(4*pow(r0,4) - 12*pow(r0,2)*pow(lambda,2) + 3*pow(lambda,4)))/(pow(E,pow(r0,2)/pow(lambda,2))*pow(lambda,4)) + 
   (352*A*pow(rr,2)*(-8*pow(r0,6) + 60*pow(r0,4)*pow(lambda,2) - 90*pow(r0,2)*pow(lambda,4) + 15*pow(lambda,6)))/(15.*pow(E,pow(r0,2)/pow(lambda,2))*pow(lambda,8)) - 
   (368*A*pow(rr,4)*(16*pow(r0,8) - 224*pow(r0,6)*pow(lambda,2) + 840*pow(r0,4)*pow(lambda,4) - 840*pow(r0,2)*pow(lambda,6) + 105*pow(lambda,8)))/
    (105.*pow(E,pow(r0,2)/pow(lambda,2))*pow(lambda,12)) + (208*A*pow(rr,6)*
      (-32*pow(r0,10) + 720*pow(r0,8)*pow(lambda,2) - 5040*pow(r0,6)*pow(lambda,4) + 12600*pow(r0,4)*pow(lambda,6) - 9450*pow(r0,2)*pow(lambda,8) + 945*pow(lambda,10)))/
    (945.*pow(E,pow(r0,2)/pow(lambda,2))*pow(lambda,16));
}
//else we continue with the numerical formulation:
else{
double F =A*(pow(E,-pow(rr - r0,2)/pow(lambda,2)) + pow(E,-pow(rr + r0,2)/pow(lambda,2)))*rr*pow(lambda,4);
double F1 = A*(pow(E,-pow(rr - r0,2)/pow(lambda,2)) + pow(E,-pow(rr + r0,2)/pow(lambda,2)))*pow(lambda,4) + 
   A*rr*((-2*(rr - r0))/(pow(E,pow(rr - r0,2)/pow(lambda,2))*pow(lambda,2)) - (2*(rr + r0))/(pow(E,pow(rr + r0,2)/pow(lambda,2))*pow(lambda,2)))*pow(lambda,4);

double F2 = A*rr*((4*pow(rr - r0,2))/(pow(E,pow(rr - r0,2)/pow(lambda,2))*pow(lambda,4)) + (4*pow(rr + r0,2))/(pow(E,pow(rr + r0,2)/pow(lambda,2))*pow(lambda,4)) - 
      2/(pow(E,pow(rr - r0,2)/pow(lambda,2))*pow(lambda,2)) - 2/(pow(E,pow(rr + r0,2)/pow(lambda,2))*pow(lambda,2)))*pow(lambda,4) + 
   2*A*((-2*(rr - r0))/(pow(E,pow(rr - r0,2)/pow(lambda,2))*pow(lambda,2)) - (2*(rr + r0))/(pow(E,pow(rr + r0,2)/pow(lambda,2))*pow(lambda,2)))*pow(lambda,4);

double F3 = A*rr*((-8*pow(rr - r0,3))/(pow(E,pow(rr - r0,2)/pow(lambda,2))*pow(lambda,6)) - (8*pow(rr + r0,3))/(pow(E,pow(rr + r0,2)/pow(lambda,2))*pow(lambda,6)) + 
      (12*(rr - r0))/(pow(E,pow(rr - r0,2)/pow(lambda,2))*pow(lambda,4)) + (12*(rr + r0))/(pow(E,pow(rr + r0,2)/pow(lambda,2))*pow(lambda,4)))*pow(lambda,4) + 
   3*A*((4*pow(rr - r0,2))/(pow(E,pow(rr - r0,2)/pow(lambda,2))*pow(lambda,4)) + (4*pow(rr + r0,2))/(pow(E,pow(rr + r0,2)/pow(lambda,2))*pow(lambda,4)) - 
      2/(pow(E,pow(rr - r0,2)/pow(lambda,2))*pow(lambda,2)) - 2/(pow(E,pow(rr + r0,2)/pow(lambda,2))*pow(lambda,2)))*pow(lambda,4);

double F4 = A*rr*((16*pow(rr - r0,4))/(pow(E,pow(rr - r0,2)/pow(lambda,2))*pow(lambda,8)) + (16*pow(rr + r0,4))/(pow(E,pow(rr + r0,2)/pow(lambda,2))*pow(lambda,8)) - 
      (48*pow(rr - r0,2))/(pow(E,pow(rr - r0,2)/pow(lambda,2))*pow(lambda,6)) - (48*pow(rr + r0,2))/(pow(E,pow(rr + r0,2)/pow(lambda,2))*pow(lambda,6)) + 
      12/(pow(E,pow(rr - r0,2)/pow(lambda,2))*pow(lambda,4)) + 12/(pow(E,pow(rr + r0,2)/pow(lambda,2))*pow(lambda,4)))*pow(lambda,4) + 
   4*A*((-8*pow(rr - r0,3))/(pow(E,pow(rr - r0,2)/pow(lambda,2))*pow(lambda,6)) - (8*pow(rr + r0,3))/(pow(E,pow(rr + r0,2)/pow(lambda,2))*pow(lambda,6)) + 
      (12*(rr - r0))/(pow(E,pow(rr - r0,2)/pow(lambda,2))*pow(lambda,4)) + (12*(rr + r0))/(pow(E,pow(rr + r0,2)/pow(lambda,2))*pow(lambda,4)))*pow(lambda,4);

A2 = 24*( -F2/(pow(rr,3))+3*F2/(pow(rr,4))-3*F/(pow(rr,5)));

B2 = 4*(-F3/(pow(rr,2))+3*F2/(pow(rr,3))- 6 * F1 /(pow(rr,4))+6*F/(pow(rr,5)));

C2 = 2*(-F4/(rr)+2*F3/(pow(rr,2))- 3* F2 /(pow(rr,3))+ 3* F1/(pow(rr,4))- 3*F/(pow(rr,5)));
}
// double A2= (48 * A * exp(-(pow((rr-r0),2)/(pow(lambda,2))))* (-2 rr* pow((rr-r0),2) - r0 * pow(lambda,2) + exp(4((rr-r0)/(pow(lambda,2)))));

// double piece1 = 48*A*(-2*rr*pow(rr + r0,2));
// double piece2 = r0*pow(lambda,2);
// double piece3 = pow(E,(4*rr*r0)/pow(lambda,2))*(-2*rr*pow(rr - r0,2));
// double piece4 = (pow(E,pow(rr + r0,2)/pow(lambda,2))*pow(rr,3));
// double B2 =(16*A*(2*pow(rr + r0,3) - 3*(rr + r0)*pow(lambda,2) + pow(E,(4*rr*r0)/pow(lambda,2))*(rr - r0)*(2*pow(rr - r0,2) - 3*pow(lambda,2))))/
//    (pow(E,pow(rr + r0,2)/pow(lambda,2))*rr*pow(lambda,2));

// double C2 = (-4*A*(pow(E,pow(rr - r0,2)/pow(lambda,2))*(8*pow(rr,3)*pow(rr + r0,4) - 8*pow(rr,2)*pow(rr + r0,2)*(4*rr + r0)*pow(lambda,2) + 6*rr*(2*pow(rr,2) - pow(r0,2))*pow(lambda,4) - 
//           3*r0*pow(lambda,6)) + pow(E,pow(rr + r0,2)/pow(lambda,2))*(8*pow(rr,3)*pow(rr - r0,4) - 8*pow(rr,2)*pow(rr - r0,2)*(4*rr - r0)*pow(lambda,2) + 
//           6*rr*(2*pow(rr,2) - pow(r0,2))*pow(lambda,4) + 3*r0*pow(lambda,6))))/(pow(E,(2*(pow(rr,2) + pow(r0,2)))/pow(lambda,2))*pow(rr,3)*pow(lambda,4));
	
// 	if (piece4 == 0) {
//     std::cout << "divide by 0!!!!" << std::endl;
// }
//     // if (isnan(A2)) printf("Error: A2 is NaN\n");
//     // if (isnan(B2)) printf("Error: B2 is NaN\n");
//     // if (isnan(C2)) printf("Error: C2 is NaN\n");
//     // if (isnan(lambda)) printf("Error: lambda is NaN\n");
//     // if (isnan(E)) printf("Error: E is NaN\n");
// 	if (isnan(piece1)) printf("Error: piece1 is NaN\n");
//     if (isnan(piece2)) printf("Error: piece2 is NaN\n");
//     if (isnan(piece3)) printf("Error: piece3 is NaN\n");
//     if (isnan(piece4)) printf("Error: piece4 is NaN\n");
//     if (isnan(E)) printf("Error: E is NaN\n");
//Defining the Spherical harmonics that we will need:
double Y20= 2- 3* sin_theta*sin_theta;
double Y20theta =-6*cos_theta*sin_theta;
double Y20_thetatheta = 3* sin_theta*sin_theta;

// //defining  metric elements in radial coordinates see equation A1
// double g_rr = (1.0+ A2*Y20 );
// double g_rth = rr*(2*B2*Y20theta);
// //double g_rph = rr* sin_theta * B2* Y20phi;
// double g_rph = 0.0;
// double g_thth = pow(rr,2)*( 1- A2*Y20/2.0 + C2*Y20_thetatheta);
// //double g_thph  = pow(rr,2)* sin_theta* (C2 * Y20_thetaphi);
// double g_thph  = 0.0;
// double g_phph  = pow(rr,2)*pow(sin_theta,2)*(1-A2*Y20/2.0 +  C2*Y20_thetatheta);



//converting into cartesian coordinates.....

//defining some things for the conversion://
// double g_bare_rr   = g_rr-1.0;
// double g_bare_thth = g_thth -rr*rr;
// double g_bare_phph = g_phph  - pow(rr,2)*pow(sin_theta,2);

double h_bare_rr   = A2* Y20;
double h_bare_rth  = (2*B2*Y20theta);
double h_bare_rph  = 0.0;
double h_bare_thth = A2*Y20/2.0 + C2*Y20_thetatheta;
double h_bare_thph = 0.0;
double h_bare_phph = A2*Y20/2.0 +  C2*Y20_thetatheta;

double g_bare_rr   = 1.0 + h_bare_rr;
double g_bare_rth  = h_bare_rth;
double g_bare_rph  = h_bare_rph;
double g_bare_thth = 1.0 - h_bare_thth;
double g_bare_thph = h_bare_thph;
double g_bare_phph = 1.0 - h_bare_phph;

double g_rr        = g_bare_rr;
double g_rth       = rr * g_bare_rth;
double g_rph       = rr * sin_theta * g_bare_rph;
double g_thth      = rr * rr * g_bare_thth;
double g_thph      = rr * rr * sin_theta * g_bare_thph;
double g_phph      = rr * rr * sin_theta * sin_theta * g_bare_phph;

if (rho > 1.0e-7) {  // everywhere but the axis ...
    drdx  = xbar / rr;
    drdy  = ybar / rr;
    drdz  = zbar / rr;

    dthdx = xbar * zbar / (rr * rr * rho);
    dthdy = ybar * zbar / (rr * rr * rho);
    dthdz = -rho / (rr * rr);

    dphdx = -ybar / (rho_sqrd);
    dphdy = xbar / (rho_sqrd);
    dphdz = 0.0;

    g_xx =
        (drdx * drdx * g_rr + dthdx * dthdx * g_thth + dphdx * dphdx * g_phph +
         2.0 * (drdx * dthdx * g_rth + drdx * dphdx * g_rph +
                dthdx * dphdx * g_thph));
    g_xy = (drdx * drdy * g_rr + dthdx * dthdy * g_thth +
            dphdx * dphdy * g_phph + (drdx * dthdy + drdy * dthdx) * g_rth +
            (drdx * dphdy + drdy * dphdx) * g_rph +
            (dthdx * dphdy + dthdy * dphdx) * g_thph);
    g_xz = (drdx * drdz * g_rr + dthdx * dthdz * g_thth +
            dphdx * dphdz * g_phph + (drdx * dthdz + drdz * dthdx) * g_rth +
            (drdx * dphdz + drdz * dphdx) * g_rph +
            (dthdx * dphdz + dthdz * dphdx) * g_thph);
    g_yy =
        (drdy * drdy * g_rr + dthdy * dthdy * g_thth + dphdy * dphdy * g_phph +
         2.0 * (drdy * dthdy * g_rth + drdy * dphdy * g_rph +
                dthdy * dphdy * g_thph));
    g_yz = (drdy * drdz * g_rr + dthdy * dthdz * g_thth +
            dphdy * dphdz * g_phph + (drdy * dthdz + drdz * dthdy) * g_rth +
            (drdy * dphdz + drdz * dphdy) * g_rph +
            (dthdy * dphdz + dthdz * dphdy) * g_thph);
    g_zz =
        (drdz * drdz * g_rr + dthdz * dthdz * g_thth + dphdz * dphdz * g_phph +
         2.0 * (drdz * dthdz * g_rth + drdz * dphdz * g_rph +
                dthdz * dphdz * g_thph));

} else if (rho < 1.0e-7 && rr > 1.0e-7) {  // on axis but not origin

    rho_axis = sqrt(rho_sqrd + epsilon_sq);

    g_xx =
        (g_bare_rr -
         (ybar * ybar + zbar * zbar) * (g_bare_rr - g_bare_thth) / (rr * rr) -
         ybar * ybar * (g_bare_thth - g_bare_phph) / (rho_axis * rho_axis) +
         2.0 / rr * (xbar / rho_axis) *
             (xbar * zbar / rr * g_bare_rth - ybar * g_bare_rph -
              zbar * (ybar / rho_axis) * g_thph));
    g_xy = (xbar * ybar / (rr * rr) * (g_bare_rr - g_bare_phph) +
            xbar * ybar / (rho_axis * rho_axis) * zbar * zbar / (rr * rr) *
                (g_bare_thth - g_bare_phph) +
            2.0 * zbar / (rr * rr) * (xbar / rho_axis) * ybar * g_bare_rth +
            (xbar / rho_axis) * xbar * g_bare_rph / rr -
            (ybar / rho_axis) * ybar * g_bare_rph / rr +
            pow((xbar / rho_axis), 2) * zbar * g_bare_thph / rr -
            pow((ybar / rho_axis), 2) * zbar * g_bare_thph / rr);
    g_xz = (xbar * zbar / (rr * rr) * (g_bare_rr - g_bare_thth) +
            (zbar * zbar * (xbar / rho_axis) * g_bare_rth -
             xbar * rho * g_bare_rth) /
                (rr * rr) -
            (zbar / rr) * (ybar / rho_axis) * g_bare_rph +
            ybar * g_bare_thph / rr);
    g_yy =
        (g_bare_rr -
         (xbar * xbar + zbar * zbar) / (rr * rr) * (g_bare_rr - g_bare_thth) -
         pow((xbar / rho_axis), 2) * (g_bare_thth - g_bare_phph) +
         2.0 / rr * (ybar / rho_axis) *
             (ybar * zbar / rr * g_bare_rth +
              (xbar / rho_axis) * zbar * g_bare_thph + xbar * g_bare_rph));
    g_yz = (ybar * zbar / (rr * rr) * (g_bare_rr - g_bare_thth) +
            (ybar / rho_axis) * zbar * zbar / (rr * rr) * g_bare_rth -
            ybar * rho / (rr * rr) * g_bare_rth +
            (zbar / rr) * (xbar / rho_axis) * g_bare_rph -
            xbar * g_bare_thph / rr);
    g_zz = (g_bare_rr - (rho * rho) / (rr * rr) * (g_bare_rr - g_bare_thth) -
            2.0 * zbar * rho * g_bare_rth / (rr * rr));

} else if (rho < 1.0e-7 && rr < 1.0e-7) {  // at the origin

    rho_axis = sqrt(rho_sqrd + epsilon_sq);
    r_origin = sqrt(r_sqrd + epsilon_sq);

    g_xx     = (g_bare_rr -
            (ybar * ybar + zbar * zbar) / (pow(r_origin, 2)) *
                (g_bare_rr - g_bare_thth) -
            ybar * ybar / (rho_axis * rho_axis) * (g_bare_thth - g_bare_phph) +
            2.0 / (r_origin) * (xbar / rho_axis) *
                (xbar * zbar / r_origin * g_bare_rth - ybar * g_bare_rph -
                 zbar * (ybar / rho_axis) * g_thph));
    g_xy     = (xbar * ybar / (pow(r_origin, 2)) * (g_bare_rr - g_bare_phph) +
            xbar * ybar / (rho_axis * rho_axis) * zbar * zbar /
                (pow(r_origin, 2)) * (g_bare_thth - g_bare_phph) +
            2.0 * zbar / (pow(r_origin, 2)) * (xbar / rho_axis) * ybar *
                g_bare_rth +
            (xbar / rho_axis) * (xbar / r_origin) * g_bare_rph -
            (ybar / rho_axis) * (ybar / r_origin) * g_bare_rph +
            pow((xbar / rho_axis), 2) * (zbar / r_origin) * g_bare_thph -
            pow((ybar / rho_axis), 2) * (zbar / r_origin) * g_bare_thph);
    g_xz     = (xbar * zbar / (pow(r_origin, 2)) * (g_bare_rr - g_bare_thth) +
            (pow((zbar / r_origin), 2) * (xbar / rho_axis) * g_bare_rth -
             (xbar / r_origin) * (rho / r_origin) * g_bare_rth) -
            (zbar / r_origin) * (ybar / rho_axis) * g_bare_rph +
            (ybar / r_origin) * g_bare_thph);
    g_yy     = (g_bare_rr -
            (pow((xbar / r_origin), 2) + pow((zbar / r_origin), 2)) *
                (g_bare_rr - g_bare_thth) -
            pow((xbar / rho_axis), 2) * (g_bare_thth - g_bare_phph) +
            2.0 * (ybar / rho_axis) *
                ((ybar / r_origin) * (zbar / r_origin) * g_bare_rth +
                 (xbar / rho_axis) * (zbar / r_origin) * g_bare_thph +
                 (xbar / r_origin) * g_bare_rph));
    g_yz = ((ybar / r_origin) * (zbar / r_origin) * (g_bare_rr - g_bare_thth) +
            (ybar / rho_axis) * pow((zbar / r_origin), 2) * g_bare_rth -
            (ybar / r_origin) * (rho / r_origin) * g_bare_rth +
            (zbar / r_origin) * (xbar / rho_axis) * g_bare_rph -
            (xbar / r_origin) * g_bare_thph);
    g_zz = (g_bare_rr - pow((rho / r_origin), 2) * (g_bare_rr - g_bare_thth) -
            2.0 * (zbar / r_origin) * (rho / r_origin) * g_bare_rth);

}

 else {
	std::cout << " something strange happened ... ";
}


// Checking to see if the determinant is negative
detg = (g_xx * g_yy * g_zz) + (g_xy * g_yz * g_xz) + (g_xz * g_xy * g_yz) -
       ((g_xz * g_yy * g_xz) + (g_yz * g_yz * g_xx) + (g_zz * g_xy * g_xy));

Chi = pow(detg, -1.0 / 3.0);
// //Check if any variable is NaN and report which one
//     if (isnan(A2)) printf("Error: A2 is NaN\n");
//     if (isnan(B2)) printf("Error: B2 is NaN\n");
//     if (isnan(C2)) printf("Error: C2 is NaN\n");
//     if (isnan(sin_theta)) printf("Error:sin is NaN\n");
//     if (isnan(cos_theta)) printf("Error: cos is NaN\n");
//     if (isnan(g_xz)) printf("Error: g_xz is NaN\n");
//     if (isnan(detg)) printf("Error: detg is NaN\n");

if (detg < 0) {
    std::cout << "The determinant = " << detg << std::endl;
    std::cout << " rr = " << rr << std::endl;
    std::cout << " x = " << xbar << std::endl;
    std::cout << " y = " << ybar << std::endl;
    std::cout << " z = " << zbar << std::endl;
	std::cout << " g_xx = " << g_xx << std::endl;
    std::cout << " g_xy = " << g_xy << std::endl;
    std::cout << " g_xz = " << g_xz << std::endl;
    std::cout << " g_yy = " << g_yy << std::endl;
    std::cout << " g_yz = " << g_yz << std::endl;
    std::cout << " g_zz = " << g_zz<< std::endl;




}

    // Other debugging options
var[VAR::U_ALPHA]  = 1.0;
var[VAR::U_BETA0]  = 0.0;
var[VAR::U_BETA1]  = 0.0;
var[VAR::U_BETA2]  = 0.0;

var[VAR::U_B0]     = 0.0;
var[VAR::U_B1]     = 0.0;
var[VAR::U_B2]     = 0.0;

var[VAR::U_GT0]    = 0.0;
var[VAR::U_GT1]    = 0.0;
var[VAR::U_GT2]    = 0.0;

var[VAR::U_SYMGT0] = Chi * g_xx;
var[VAR::U_SYMGT1] = Chi * g_xy;
var[VAR::U_SYMGT2] = Chi * g_xz;
var[VAR::U_SYMGT3] = Chi * g_yy;
var[VAR::U_SYMGT4] = Chi * g_yz;
var[VAR::U_SYMGT5] = Chi * g_zz;
var[VAR::U_SYMAT0] = 0.0;
var[VAR::U_SYMAT1] = 0.0;
var[VAR::U_SYMAT2] = 0.0;
var[VAR::U_SYMAT3] = 0.0;
var[VAR::U_SYMAT4] = 0.0;
var[VAR::U_SYMAT5] = 0.0;

var[VAR::U_CHI]    = Chi;
var[VAR::U_K]      = 0.0;
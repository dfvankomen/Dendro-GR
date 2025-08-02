#include "Linear_Teukolsky.h"

namespace bssn {
void LinearTeuk(const double xx1, const double yy1, const double zz1,
                const double t, double *var, bool varsAreGrid) {
    double xx, yy, zz;
    if (varsAreGrid) {
        xx = GRIDX_TO_X(xx1);
        yy = GRIDY_TO_Y(yy1);
        zz = GRIDZ_TO_Z(zz1);
    } else {
        xx = xx1;
        yy = yy1;
        zz = zz1;
    }
//This is an initial data file design to match the one used by Baumgarte: https://arxiv.org/abs/2303.05530
//reading in data from par file...
double A = bssn::TEUK_AMP;
double lambda = bssn::TEUK_WIDTH; //half width really; width is twice this
double r0 = bssn::TEUK_R_0;
//These are given in a mathematica script and converted into C code:
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
int teukolsky_wave_gauge;
int multipole_l;
int multipole_m;
int multipole_parity;
int id_dynamics_type;
int id_functional_form;
int kk;
double sin_theta;
double cos_theta;
double sin_2theta;
double cos_2theta;
double sin_3theta;
double cos_3theta;
double sin_4theta;
double cos_4theta;
double sin_phi;
double cos_phi;
double sin_2phi;
double cos_2phi;
double sin_3phi;
double cos_3phi;
double sin_4phi;
double cos_4phi;
//double rho_sqrd;
//double rho;
//double rho_axis;
//double rr;
//double r_origin;
//double r_sqrd;
double a_0;
double b_0;
double r_0;
double r_0_sqrd;
double width_teuk;
double aa;
double bb;
double amp_teuk;
double f0;
double df0;
double ddf0;
double d3f0;
double d4f0;
double d5f0;
double d6f0;
double f1;
double f2;
double df1;
double df2;
double ddf1;
double ddf2;
double d3f1;
double d3f2;
double d4f1;
double d4f2;
double d5f1;
double d5f2;
double d6f1;
double d6f2;
double A_lm;
double B_lm;
double C_lm;
double K_lm;
double L_lm;
// double drdx;
// double drdy;
// double drdz;
// double dthdx;
// double dthdy;
// double dthdz;
// double dphdx;
// double dphdy;
// double dphdz;
// double g_xx;
// double g_xy;
// double g_xz;
// double g_yy;
// double g_yz;
// double g_zz;
//double epsilon_sq;
double Ylm;
double Ylm_theta;
double Ylm_phi;
double Ylm_thth;
double Ylm_thphi;
double Y20;
double Y21;
double Y2n1;
double Y22;
double Y2n2;
double Y20_theta;
double Y21_theta;
double Y2n1_theta;
double Y22_theta;
double Y2n2_theta;
double Y20_phi;
double Y21_phi;
double Y2n1_phi;
double Y22_phi;
double Y2n2_phi;
double Y20_thth;
double Y21_thth;
double Y2n1_thth;
double Y22_thth;
double Y2n2_thth;
double Y20_thphi;
double Y21_thphi;
double Y2n1_thphi;
double Y22_thphi;
double Y2n2_thphi;
double Y30;
double Y31;
double Y3n1;
double Y32;
double Y3n2;
double Y33;
double Y3n3;
double Y30_theta;
double Y31_theta;
double Y3n1_theta;
double Y32_theta;
double Y3n2_theta;
double Y33_theta;
double Y3n3_theta;
double Y30_phi;
double Y31_phi;
double Y3n1_phi;
double Y32_phi;
double Y3n2_phi;
double Y33_phi;
double Y3n3_phi;
double Y30_thth;
double Y31_thth;
double Y3n1_thth;
double Y32_thth;
double Y3n2_thth;
double Y33_thth;
double Y3n3_thth;
double Y30_thphi;
double Y31_thphi;
double Y3n1_thphi;
double Y32_thphi;
double Y3n2_thphi;
double Y33_thphi;
double Y3n3_thphi;
double Y40;
double Y41;
double Y4n1;
double Y42;
double Y4n2;
double Y43;
double Y4n3;
double Y44;
double Y4n4;
double Y40_theta;
double Y41_theta;
double Y4n1_theta;
double Y42_theta;
double Y4n2_theta;
double Y43_theta;
double Y4n3_theta;
double Y44_theta;
double Y4n4_theta;
double Y40_phi;
double Y41_phi;
double Y4n1_phi;
double Y42_phi;
double Y4n2_phi;
double Y43_phi;
double Y4n3_phi;
double Y44_phi;
double Y4n4_phi;
double Y40_thth;
double Y41_thth;
double Y4n1_thth;
double Y42_thth;
double Y4n2_thth;
double Y43_thth;
double Y4n3_thth;
double Y44_thth;
double Y4n4_thth;
double Y40_thphi;
double Y41_thphi;
double Y4n1_thphi;
double Y42_thphi;
double Y4n2_thphi;
double Y43_thphi;
double Y4n3_thphi;
double Y44_thphi;
double Y4n4_thphi;
double h_bare_rr;
double h_bare_rth;
double h_bare_rph;
double h_bare_thth;
double h_bare_thph;
double h_bare_phph;
double g_bare_rr;
double g_bare_rth;
double g_bare_rph;
double g_bare_thth;
double g_bare_thph;
double g_bare_phph;
double g_rr;
double g_rth;
double g_rph;
double g_thth;
double g_thph;
double g_phph;
//double detg;
//double pi;
double chi;
double fr0;

// std::cout << "starting the loop in initial data ... /n";
//double offset = 0.01627605 / 7.0;
double xbar   = xx;
double ybar   = yy;
double zbar   = zz;

double rho_sqrd      = xbar * xbar + ybar * ybar;
double rho           = sqrt(rho_sqrd);

double r_sqrd        = rho_sqrd + zbar * zbar;
double rr            = sqrt(r_sqrd);

double epsilon_sq    = 1.0e-5;  // maybe we can read in grid info later?
double r_origin      = sqrt(r_sqrd + epsilon_sq);
if(rr<1e-6){
    rr=r_origin;
}
// Away from the axis of symmetry
if (abs(xbar) > 1.0e-7 || abs(ybar) > 1.0e-7) {
    sin_theta  = rho / rr;
    cos_theta  = zbar / rr;
    sin_2theta = 2.0 * sin_theta * cos_theta;
    cos_2theta = pow(cos_theta, 2) - pow(sin_theta, 2);
    sin_3theta = (3.0 * sin_theta) - (4.0 * pow(sin_theta, 3));
    cos_3theta = (4.0 * pow(cos_theta, 3)) - (3.0 * cos_theta);
    sin_4theta = cos_theta * (4.0 * sin_theta - 8.0 * pow(sin_theta, 3));
    cos_4theta = 8.0 * pow(cos_theta, 4) - 8.0 * pow(cos_theta, 2) + 1.0;
    sin_phi    = ybar / rho;
    cos_phi    = xbar / rho;
    sin_2phi   = 2.0 * sin_phi * cos_phi;
    cos_2phi   = pow(cos_phi, 2) - pow(sin_phi, 2);
    sin_3phi   = (3.0 * sin_phi) - (4.0 * pow(sin_phi, 3));
    cos_3phi   = (4.0 * pow(cos_phi, 3)) - (3.0 * cos_phi);
    sin_4phi   = cos_phi * (4.0 * sin_phi - 8.0 * pow(sin_phi, 3));
    cos_4phi   = 8.0 * pow(cos_phi, 4) - 8.0 * pow(cos_phi, 2) + 1.0;
}
// On the axis of symmetry but away from the origin
else if ((abs(xbar) < 1.0e-7 && abs(ybar) < 1.0e-7) && abs(zbar) > 1.0e-7) {
    sin_theta  = rho / rr;
    cos_theta  = zbar / rr;
    sin_2theta = 2.0 * sin_theta * cos_theta;
    cos_2theta = pow(cos_theta, 2) - pow(sin_theta, 2);
    sin_3theta = (3.0 * sin_theta) - (4.0 * pow(sin_theta, 3));
    cos_3theta = (4.0 * pow(cos_theta, 3)) - (3.0 * cos_theta);
    sin_4theta = cos_theta * (4.0 * sin_theta - 8.0 * pow(sin_theta, 3));
    cos_4theta = 8.0 * pow(cos_theta, 4) - 8.0 * pow(cos_theta, 2) + 1.0;
    sin_phi    = 0.0;
    cos_phi    = 0.0;
    sin_2phi   = 2.0 * sin_phi * cos_phi;
    cos_2phi   = pow(cos_phi, 2) - pow(sin_phi, 2);
    sin_3phi   = (3.0 * sin_phi) - (4.0 * pow(sin_phi, 3));
    cos_3phi   = (4.0 * pow(cos_phi, 3)) - (3.0 * cos_phi);
    sin_4phi   = cos_phi * (4.0 * sin_phi - 8.0 * pow(sin_phi, 3));
    cos_4phi   = 8.0 * pow(cos_phi, 4) - 8.0 * pow(cos_phi, 2) + 1.0;
}
// At the origin
else if ((abs(xbar) < 1.0e-7 && abs(ybar) < 1.0e-7) && abs(zbar) < 1.0e-7) {
    sin_theta  = 0.0;
    cos_theta  = 0.0;
    sin_2theta = 2.0 * sin_theta * cos_theta;
    cos_2theta = pow(cos_theta, 2) - pow(sin_theta, 2);
    sin_3theta = (3.0 * sin_theta) - (4.0 * pow(sin_theta, 3));
    cos_3theta = (4.0 * pow(cos_theta, 3)) - (3.0 * cos_theta);
    sin_4theta = cos_theta * (4.0 * sin_theta - 8.0 * pow(sin_theta, 3));
    cos_4theta = 8.0 * pow(cos_theta, 4) - 8.0 * pow(cos_theta, 2) + 1.0;
    sin_phi    = 0.0;
    cos_phi    = 0.0;
    sin_2phi   = 2.0 * sin_phi * cos_phi;
    cos_2phi   = pow(cos_phi, 2) - pow(sin_phi, 2);
    sin_3phi   = (3.0 * sin_phi) - (4.0 * pow(sin_phi, 3));
    cos_3phi   = (4.0 * pow(cos_phi, 3)) - (3.0 * cos_phi);
    sin_4phi   = cos_phi * (4.0 * sin_phi - 8.0 * pow(sin_phi, 3));
    cos_4phi   = 8.0 * pow(cos_phi, 4) - 8.0 * pow(cos_phi, 2) + 1.0;
} else {
    std::cout << " problem in initial data /n";
    std::cout << " x = " << xbar;
    std::cout << " y = " << ybar;
    std::cout << " z = " << zbar;
    std::cout << " r = " << rr;
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
//we will here make use of a taylor expanded expression of these functions as to avoid singularities arising from imperfect cancelations at the origin  We do this out to order 6. 
//At some point if I or some other unfortunate soul has to work on this I will note here that if r0 or t get too large then this will fail at the origin as we will have exp(-big number) so this will then fail. There is probably an intelligent wway to drop these terms but I haven't taken the time to figure it out. 
if (rr < 0.75 ) {
A2 = (72*A*(-1 + lambda)*(-2*pow(t,2) - 2*r0*t*lambda + pow(lambda,4) + pow(E,(4*r0*t)/pow(lambda,3))*(-2*pow(t,2) + 2*r0*t*lambda + pow(lambda,4))))/(pow(E,pow(t + r0*lambda,2)/pow(lambda,4))*pow(rr,4)) + 
   (24*A*(-1 + lambda)*(-1 + 2*lambda)*(pow(E,pow(t - r0*lambda,2)/pow(lambda,4))*(4*pow(t,4) + 12*r0*pow(t,3)*lambda - 6*pow(r0,2)*pow(lambda,6) + 3*pow(lambda,8) + 12*pow(t,2)*(r0 - lambda)*pow(lambda,2)*(r0 + lambda) + 
           2*t*(2*pow(r0,3)*pow(lambda,3) - 9*r0*pow(lambda,5))) + pow(E,pow(t + r0*lambda,2)/pow(lambda,4))*
         (4*pow(t,4) - 12*r0*pow(t,3)*lambda - 6*pow(r0,2)*pow(lambda,6) + 3*pow(lambda,8) + 12*pow(t,2)*(r0 - lambda)*pow(lambda,2)*(r0 + lambda) + t*(-4*pow(r0,3)*pow(lambda,3) + 18*r0*pow(lambda,5)))))/
    (pow(E,(2*(pow(t,2) + pow(r0,2)*pow(lambda,2)))/pow(lambda,4))*pow(rr,2)*pow(lambda,8)) - (A*pow(rr,4)*(1 + 3*lambda*(-3 + 8*lambda))*
      (pow(E,pow(t - r0*lambda,2)/pow(lambda,4))*(-32*pow(t,10) - 288*r0*pow(t,9)*lambda + 144*pow(t,8)*pow(lambda,2)*(-8*pow(r0,2) + 5*pow(lambda,2)) + 192*r0*pow(t,7)*pow(lambda,3)*(-14*pow(r0,2) + 27*pow(lambda,2)) - 
           144*r0*pow(t,3)*pow(lambda,7)*(8*pow(r0,6) - 140*pow(r0,4)*pow(lambda,2) + 490*pow(r0,2)*pow(lambda,4) - 315*pow(lambda,6)) + 
           168*pow(t,4)*pow(lambda,6)*(-16*pow(r0,6) + 180*pow(r0,4)*pow(lambda,2) - 360*pow(r0,2)*pow(lambda,4) + 75*pow(lambda,6)) - 1008*pow(t,6)*(4*pow(r0,4)*pow(lambda,4) - 16*pow(r0,2)*pow(lambda,6) + 5*pow(lambda,8)) + 
           9*pow(lambda,12)*(16*pow(r0,8) - 224*pow(r0,6)*pow(lambda,2) + 840*pow(r0,4)*pow(lambda,4) - 840*pow(r0,2)*pow(lambda,6) + 105*pow(lambda,8)) - 
           18*pow(t,2)*pow(lambda,8)*(16*pow(r0,8) - 448*pow(r0,6)*pow(lambda,2) + 2520*pow(r0,4)*pow(lambda,4) - 3360*pow(r0,2)*pow(lambda,6) + 525*pow(lambda,8)) - 
           2*r0*t*pow(lambda,9)*(16*pow(r0,8) - 864*pow(r0,6)*pow(lambda,2) + 7560*pow(r0,4)*pow(lambda,4) - 17640*pow(r0,2)*pow(lambda,6) + 8505*pow(lambda,8)) - 
           1008*pow(t,5)*(4*pow(r0,5)*pow(lambda,5) - 28*pow(r0,3)*pow(lambda,7) + 27*r0*pow(lambda,9))) + 
        pow(E,pow(t + r0*lambda,2)/pow(lambda,4))*(-32*pow(t,10) + 288*r0*pow(t,9)*lambda + 144*pow(t,8)*pow(lambda,2)*(-8*pow(r0,2) + 5*pow(lambda,2)) + 192*pow(t,7)*(14*pow(r0,3)*pow(lambda,3) - 27*r0*pow(lambda,5)) + 
           144*r0*pow(t,3)*pow(lambda,7)*(8*pow(r0,6) - 140*pow(r0,4)*pow(lambda,2) + 490*pow(r0,2)*pow(lambda,4) - 315*pow(lambda,6)) + 
           168*pow(t,4)*pow(lambda,6)*(-16*pow(r0,6) + 180*pow(r0,4)*pow(lambda,2) - 360*pow(r0,2)*pow(lambda,4) + 75*pow(lambda,6)) - 1008*pow(t,6)*(4*pow(r0,4)*pow(lambda,4) - 16*pow(r0,2)*pow(lambda,6) + 5*pow(lambda,8)) + 
           9*pow(lambda,12)*(16*pow(r0,8) - 224*pow(r0,6)*pow(lambda,2) + 840*pow(r0,4)*pow(lambda,4) - 840*pow(r0,2)*pow(lambda,6) + 105*pow(lambda,8)) - 
           18*pow(t,2)*pow(lambda,8)*(16*pow(r0,8) - 448*pow(r0,6)*pow(lambda,2) + 2520*pow(r0,4)*pow(lambda,4) - 3360*pow(r0,2)*pow(lambda,6) + 525*pow(lambda,8)) + 
           2*r0*t*pow(lambda,9)*(16*pow(r0,8) - 864*pow(r0,6)*pow(lambda,2) + 7560*pow(r0,4)*pow(lambda,4) - 17640*pow(r0,2)*pow(lambda,6) + 8505*pow(lambda,8)) + 
           1008*pow(t,5)*(4*pow(r0,5)*pow(lambda,5) - 28*pow(r0,3)*pow(lambda,7) + 27*r0*pow(lambda,9)))))/(315.*pow(E,(2*(pow(t,2) + pow(r0,2)*pow(lambda,2)))/pow(lambda,4))*pow(lambda,32)) - 
   (4*A*(3 + 5*lambda*(-3 + 4*lambda))*(pow(E,pow(t - r0*lambda,2)/pow(lambda,4))*(-8*pow(t,6) - 40*r0*pow(t,5)*lambda + 40*r0*pow(t,3)*pow(lambda,3)*(-2*pow(r0,2) + 5*pow(lambda,2)) + 
           5*pow(lambda,8)*(4*pow(r0,4) - 12*pow(r0,2)*pow(lambda,2) + 3*pow(lambda,4)) + pow(t,4)*(-80*pow(r0,2)*pow(lambda,2) + 60*pow(lambda,4)) - 10*pow(t,2)*(4*pow(r0,4)*pow(lambda,4) - 24*pow(r0,2)*pow(lambda,6) + 9*pow(lambda,8)) - 
           2*t*(4*pow(r0,5)*pow(lambda,5) - 60*pow(r0,3)*pow(lambda,7) + 75*r0*pow(lambda,9))) + pow(E,pow(t + r0*lambda,2)/pow(lambda,4))*
         (-8*pow(t,6) + 40*r0*pow(t,5)*lambda + 5*pow(lambda,8)*(4*pow(r0,4) - 12*pow(r0,2)*pow(lambda,2) + 3*pow(lambda,4)) + pow(t,4)*(-80*pow(r0,2)*pow(lambda,2) + 60*pow(lambda,4)) + 
           40*pow(t,3)*(2*pow(r0,3)*pow(lambda,3) - 5*r0*pow(lambda,5)) - 10*pow(t,2)*(4*pow(r0,4)*pow(lambda,4) - 24*pow(r0,2)*pow(lambda,6) + 9*pow(lambda,8)) + 2*t*(4*pow(r0,5)*pow(lambda,5) - 60*pow(r0,3)*pow(lambda,7) + 75*r0*pow(lambda,9)))
        ))/(5.*pow(E,(2*(pow(t,2) + pow(r0,2)*pow(lambda,2)))/pow(lambda,4))*pow(lambda,16)) + (4*A*pow(rr,2)*(1 + 7*lambda*(-1 + 2*lambda))*
      (pow(E,pow(t + r0*lambda,2)/pow(lambda,4))*(16*pow(t,8) - 112*r0*pow(t,7)*lambda + 56*r0*pow(t,5)*pow(lambda,3)*(-10*pow(r0,2) + 21*pow(lambda,2)) + 112*pow(t,6)*(3*pow(r0,2)*pow(lambda,2) - 2*pow(lambda,4)) + 
           28*pow(t,2)*pow(lambda,6)*(4*pow(r0,6) - 60*pow(r0,4)*pow(lambda,2) + 135*pow(r0,2)*pow(lambda,4) - 30*pow(lambda,6)) + 7*pow(lambda,10)*(-8*pow(r0,6) + 60*pow(r0,4)*pow(lambda,2) - 90*pow(r0,2)*pow(lambda,4) + 15*pow(lambda,6)) + 
           2*r0*t*pow(lambda,7)*(-8*pow(r0,6) + 252*pow(r0,4)*pow(lambda,2) - 1050*pow(r0,2)*pow(lambda,4) + 735*pow(lambda,6)) + 280*pow(t,4)*(2*pow(r0,4)*pow(lambda,4) - 9*pow(r0,2)*pow(lambda,6) + 3*pow(lambda,8)) - 
           28*pow(t,3)*(12*pow(r0,5)*pow(lambda,5) - 100*pow(r0,3)*pow(lambda,7) + 105*r0*pow(lambda,9))) + 
        pow(E,pow(t - r0*lambda,2)/pow(lambda,4))*(16*pow(t,8) + 112*r0*pow(t,7)*lambda + 112*pow(t,6)*(3*pow(r0,2)*pow(lambda,2) - 2*pow(lambda,4)) + 56*pow(t,5)*(10*pow(r0,3)*pow(lambda,3) - 21*r0*pow(lambda,5)) + 
           2*r0*t*pow(lambda,7)*(8*pow(r0,6) - 252*pow(r0,4)*pow(lambda,2) + 1050*pow(r0,2)*pow(lambda,4) - 735*pow(lambda,6)) + 
           28*pow(t,2)*pow(lambda,6)*(4*pow(r0,6) - 60*pow(r0,4)*pow(lambda,2) + 135*pow(r0,2)*pow(lambda,4) - 30*pow(lambda,6)) + 7*pow(lambda,10)*(-8*pow(r0,6) + 60*pow(r0,4)*pow(lambda,2) - 90*pow(r0,2)*pow(lambda,4) + 15*pow(lambda,6)) + 
           280*pow(t,4)*(2*pow(r0,4)*pow(lambda,4) - 9*pow(r0,2)*pow(lambda,6) + 3*pow(lambda,8)) + 28*pow(t,3)*(12*pow(r0,5)*pow(lambda,5) - 100*pow(r0,3)*pow(lambda,7) + 105*r0*pow(lambda,9)))))/
    (35.*pow(E,(2*(pow(t,2) + pow(r0,2)*pow(lambda,2)))/pow(lambda,4))*pow(lambda,24)) + (A*pow(rr,6)*(3 + 11*lambda*(-3 + 10*lambda))*
      (pow(E,pow(t - r0*lambda,2)/pow(lambda,4))*(64*pow(t,12) + 704*r0*pow(t,11)*lambda + 704*pow(t,10)*(5*pow(r0,2)*pow(lambda,2) - 3*pow(lambda,4)) + 1760*pow(t,9)*(6*pow(r0,3)*pow(lambda,3) - 11*r0*pow(lambda,5)) + 
           528*r0*pow(t,5)*pow(lambda,7)*(40*pow(r0,6) - 588*pow(r0,4)*pow(lambda,2) + 1890*pow(r0,2)*pow(lambda,4) - 1155*pow(lambda,6)) + 
           7392*pow(t,6)*pow(lambda,6)*(4*pow(r0,6) - 40*pow(r0,4)*pow(lambda,2) + 75*pow(r0,2)*pow(lambda,4) - 15*pow(lambda,6)) + 2640*pow(t,8)*(8*pow(r0,4)*pow(lambda,4) - 30*pow(r0,2)*pow(lambda,6) + 9*pow(lambda,8)) + 
           660*pow(t,4)*pow(lambda,8)*(16*pow(r0,8) - 336*pow(r0,6)*pow(lambda,2) + 1680*pow(r0,4)*pow(lambda,4) - 2100*pow(r0,2)*pow(lambda,6) + 315*pow(lambda,8)) + 
           220*r0*pow(t,3)*pow(lambda,9)*(16*pow(r0,8) - 480*pow(r0,6)*pow(lambda,2) + 3528*pow(r0,4)*pow(lambda,4) - 7560*pow(r0,2)*pow(lambda,6) + 3465*pow(lambda,8)) + 
           1056*pow(t,7)*(28*pow(r0,5)*pow(lambda,5) - 180*pow(r0,3)*pow(lambda,7) + 165*r0*pow(lambda,9)) + 
           2*r0*t*pow(lambda,11)*(32*pow(r0,10) - 2640*pow(r0,8)*pow(lambda,2) + 39600*pow(r0,6)*pow(lambda,4) - 194040*pow(r0,4)*pow(lambda,6) + 311850*pow(r0,2)*pow(lambda,8) - 114345*pow(lambda,10)) + 
           44*pow(t,2)*pow(lambda,10)*(16*pow(r0,10) - 720*pow(r0,8)*pow(lambda,2) + 7560*pow(r0,6)*pow(lambda,4) - 25200*pow(r0,4)*pow(lambda,6) + 23625*pow(r0,2)*pow(lambda,8) - 2835*pow(lambda,10)) + 
           11*pow(lambda,14)*(-32*pow(r0,10) + 720*pow(r0,8)*pow(lambda,2) - 5040*pow(r0,6)*pow(lambda,4) + 12600*pow(r0,4)*pow(lambda,6) - 9450*pow(r0,2)*pow(lambda,8) + 945*pow(lambda,10))) + 
        pow(E,pow(t + r0*lambda,2)/pow(lambda,4))*(64*pow(t,12) - 704*r0*pow(t,11)*lambda + 704*pow(t,10)*(5*pow(r0,2)*pow(lambda,2) - 3*pow(lambda,4)) - 1760*pow(t,9)*(6*pow(r0,3)*pow(lambda,3) - 11*r0*pow(lambda,5)) - 
           528*r0*pow(t,5)*pow(lambda,7)*(40*pow(r0,6) - 588*pow(r0,4)*pow(lambda,2) + 1890*pow(r0,2)*pow(lambda,4) - 1155*pow(lambda,6)) + 
           7392*pow(t,6)*pow(lambda,6)*(4*pow(r0,6) - 40*pow(r0,4)*pow(lambda,2) + 75*pow(r0,2)*pow(lambda,4) - 15*pow(lambda,6)) + 2640*pow(t,8)*(8*pow(r0,4)*pow(lambda,4) - 30*pow(r0,2)*pow(lambda,6) + 9*pow(lambda,8)) + 
           660*pow(t,4)*pow(lambda,8)*(16*pow(r0,8) - 336*pow(r0,6)*pow(lambda,2) + 1680*pow(r0,4)*pow(lambda,4) - 2100*pow(r0,2)*pow(lambda,6) + 315*pow(lambda,8)) - 
           220*r0*pow(t,3)*pow(lambda,9)*(16*pow(r0,8) - 480*pow(r0,6)*pow(lambda,2) + 3528*pow(r0,4)*pow(lambda,4) - 7560*pow(r0,2)*pow(lambda,6) + 3465*pow(lambda,8)) - 
           1056*pow(t,7)*(28*pow(r0,5)*pow(lambda,5) - 180*pow(r0,3)*pow(lambda,7) + 165*r0*pow(lambda,9)) + 
           44*pow(t,2)*pow(lambda,10)*(16*pow(r0,10) - 720*pow(r0,8)*pow(lambda,2) + 7560*pow(r0,6)*pow(lambda,4) - 25200*pow(r0,4)*pow(lambda,6) + 23625*pow(r0,2)*pow(lambda,8) - 2835*pow(lambda,10)) + 
           11*pow(lambda,14)*(-32*pow(r0,10) + 720*pow(r0,8)*pow(lambda,2) - 5040*pow(r0,6)*pow(lambda,4) + 12600*pow(r0,4)*pow(lambda,6) - 9450*pow(r0,2)*pow(lambda,8) + 945*pow(lambda,10)) + 
           2*r0*t*pow(lambda,11)*(-32*pow(r0,10) + 2640*pow(r0,8)*pow(lambda,2) - 39600*pow(r0,6)*pow(lambda,4) + 194040*pow(r0,4)*pow(lambda,6) - 311850*pow(r0,2)*pow(lambda,8) + 114345*pow(lambda,10)))))/
    (51975.*pow(E,(2*(pow(t,2) + pow(r0,2)*pow(lambda,2)))/pow(lambda,4))*pow(lambda,40));
B2= (24*A*(-1 + lambda)*(2*pow(t,2) + 2*r0*t*lambda - pow(lambda,4) - pow(E,(4*r0*t)/pow(lambda,3))*(-2*pow(t,2) + 2*r0*t*lambda + pow(lambda,4))))/(pow(E,pow(t + r0*lambda,2)/pow(lambda,4))*pow(rr,4)) + 
   (8*A*pow(-1 + lambda,3)*(pow(E,pow(t - r0*lambda,2)/pow(lambda,4))*(4*pow(t,4) + 12*r0*pow(t,3)*lambda - 6*pow(r0,2)*pow(lambda,6) + 3*pow(lambda,8) + 12*pow(t,2)*(r0 - lambda)*pow(lambda,2)*(r0 + lambda) + 
           2*t*(2*pow(r0,3)*pow(lambda,3) - 9*r0*pow(lambda,5))) + pow(E,pow(t + r0*lambda,2)/pow(lambda,4))*
         (4*pow(t,4) - 12*r0*pow(t,3)*lambda - 6*pow(r0,2)*pow(lambda,6) + 3*pow(lambda,8) + 12*pow(t,2)*(r0 - lambda)*pow(lambda,2)*(r0 + lambda) + t*(-4*pow(r0,3)*pow(lambda,3) + 18*r0*pow(lambda,5)))))/
    (pow(E,(2*(pow(t,2) + pow(r0,2)*pow(lambda,2)))/pow(lambda,4))*pow(rr,2)*pow(lambda,8)) - (A*pow(rr,4)*(-1 + 3*lambda*(3 + 4*lambda*(-3 + 7*lambda)))*
      (pow(E,pow(t - r0*lambda,2)/pow(lambda,4))*(-32*pow(t,10) - 288*r0*pow(t,9)*lambda + 144*pow(t,8)*pow(lambda,2)*(-8*pow(r0,2) + 5*pow(lambda,2)) + 192*r0*pow(t,7)*pow(lambda,3)*(-14*pow(r0,2) + 27*pow(lambda,2)) - 
           144*r0*pow(t,3)*pow(lambda,7)*(8*pow(r0,6) - 140*pow(r0,4)*pow(lambda,2) + 490*pow(r0,2)*pow(lambda,4) - 315*pow(lambda,6)) + 
           168*pow(t,4)*pow(lambda,6)*(-16*pow(r0,6) + 180*pow(r0,4)*pow(lambda,2) - 360*pow(r0,2)*pow(lambda,4) + 75*pow(lambda,6)) - 1008*pow(t,6)*(4*pow(r0,4)*pow(lambda,4) - 16*pow(r0,2)*pow(lambda,6) + 5*pow(lambda,8)) + 
           9*pow(lambda,12)*(16*pow(r0,8) - 224*pow(r0,6)*pow(lambda,2) + 840*pow(r0,4)*pow(lambda,4) - 840*pow(r0,2)*pow(lambda,6) + 105*pow(lambda,8)) - 
           18*pow(t,2)*pow(lambda,8)*(16*pow(r0,8) - 448*pow(r0,6)*pow(lambda,2) + 2520*pow(r0,4)*pow(lambda,4) - 3360*pow(r0,2)*pow(lambda,6) + 525*pow(lambda,8)) - 
           2*r0*t*pow(lambda,9)*(16*pow(r0,8) - 864*pow(r0,6)*pow(lambda,2) + 7560*pow(r0,4)*pow(lambda,4) - 17640*pow(r0,2)*pow(lambda,6) + 8505*pow(lambda,8)) - 
           1008*pow(t,5)*(4*pow(r0,5)*pow(lambda,5) - 28*pow(r0,3)*pow(lambda,7) + 27*r0*pow(lambda,9))) + 
        pow(E,pow(t + r0*lambda,2)/pow(lambda,4))*(-32*pow(t,10) + 288*r0*pow(t,9)*lambda + 144*pow(t,8)*pow(lambda,2)*(-8*pow(r0,2) + 5*pow(lambda,2)) + 192*pow(t,7)*(14*pow(r0,3)*pow(lambda,3) - 27*r0*pow(lambda,5)) + 
           144*r0*pow(t,3)*pow(lambda,7)*(8*pow(r0,6) - 140*pow(r0,4)*pow(lambda,2) + 490*pow(r0,2)*pow(lambda,4) - 315*pow(lambda,6)) + 
           168*pow(t,4)*pow(lambda,6)*(-16*pow(r0,6) + 180*pow(r0,4)*pow(lambda,2) - 360*pow(r0,2)*pow(lambda,4) + 75*pow(lambda,6)) - 1008*pow(t,6)*(4*pow(r0,4)*pow(lambda,4) - 16*pow(r0,2)*pow(lambda,6) + 5*pow(lambda,8)) + 
           9*pow(lambda,12)*(16*pow(r0,8) - 224*pow(r0,6)*pow(lambda,2) + 840*pow(r0,4)*pow(lambda,4) - 840*pow(r0,2)*pow(lambda,6) + 105*pow(lambda,8)) - 
           18*pow(t,2)*pow(lambda,8)*(16*pow(r0,8) - 448*pow(r0,6)*pow(lambda,2) + 2520*pow(r0,4)*pow(lambda,4) - 3360*pow(r0,2)*pow(lambda,6) + 525*pow(lambda,8)) + 
           2*r0*t*pow(lambda,9)*(16*pow(r0,8) - 864*pow(r0,6)*pow(lambda,2) + 7560*pow(r0,4)*pow(lambda,4) - 17640*pow(r0,2)*pow(lambda,6) + 8505*pow(lambda,8)) + 
           1008*pow(t,5)*(4*pow(r0,5)*pow(lambda,5) - 28*pow(r0,3)*pow(lambda,7) + 27*r0*pow(lambda,9)))))/(945.*pow(E,(2*(pow(t,2) + pow(r0,2)*pow(lambda,2)))/pow(lambda,4))*pow(lambda,32)) - 
   (4*A*(-1 + 5*lambda*(1 + 2*(-1 + lambda)*lambda))*(pow(E,pow(t - r0*lambda,2)/pow(lambda,4))*(-8*pow(t,6) - 40*r0*pow(t,5)*lambda + 40*r0*pow(t,3)*pow(lambda,3)*(-2*pow(r0,2) + 5*pow(lambda,2)) + 
           5*pow(lambda,8)*(4*pow(r0,4) - 12*pow(r0,2)*pow(lambda,2) + 3*pow(lambda,4)) + pow(t,4)*(-80*pow(r0,2)*pow(lambda,2) + 60*pow(lambda,4)) - 10*pow(t,2)*(4*pow(r0,4)*pow(lambda,4) - 24*pow(r0,2)*pow(lambda,6) + 9*pow(lambda,8)) - 
           2*t*(4*pow(r0,5)*pow(lambda,5) - 60*pow(r0,3)*pow(lambda,7) + 75*r0*pow(lambda,9))) + pow(E,pow(t + r0*lambda,2)/pow(lambda,4))*
         (-8*pow(t,6) + 40*r0*pow(t,5)*lambda + 5*pow(lambda,8)*(4*pow(r0,4) - 12*pow(r0,2)*pow(lambda,2) + 3*pow(lambda,4)) + pow(t,4)*(-80*pow(r0,2)*pow(lambda,2) + 60*pow(lambda,4)) + 
           40*pow(t,3)*(2*pow(r0,3)*pow(lambda,3) - 5*r0*pow(lambda,5)) - 10*pow(t,2)*(4*pow(r0,4)*pow(lambda,4) - 24*pow(r0,2)*pow(lambda,6) + 9*pow(lambda,8)) + 2*t*(4*pow(r0,5)*pow(lambda,5) - 60*pow(r0,3)*pow(lambda,7) + 75*r0*pow(lambda,9)))
        ))/(5.*pow(E,(2*(pow(t,2) + pow(r0,2)*pow(lambda,2)))/pow(lambda,4))*pow(lambda,16)) + (4*A*pow(rr,2)*(-1 + 7*lambda*(1 + lambda*(-3 + 5*lambda)))*
      (pow(E,pow(t + r0*lambda,2)/pow(lambda,4))*(16*pow(t,8) - 112*r0*pow(t,7)*lambda + 56*r0*pow(t,5)*pow(lambda,3)*(-10*pow(r0,2) + 21*pow(lambda,2)) + 112*pow(t,6)*(3*pow(r0,2)*pow(lambda,2) - 2*pow(lambda,4)) + 
           28*pow(t,2)*pow(lambda,6)*(4*pow(r0,6) - 60*pow(r0,4)*pow(lambda,2) + 135*pow(r0,2)*pow(lambda,4) - 30*pow(lambda,6)) + 7*pow(lambda,10)*(-8*pow(r0,6) + 60*pow(r0,4)*pow(lambda,2) - 90*pow(r0,2)*pow(lambda,4) + 15*pow(lambda,6)) + 
           2*r0*t*pow(lambda,7)*(-8*pow(r0,6) + 252*pow(r0,4)*pow(lambda,2) - 1050*pow(r0,2)*pow(lambda,4) + 735*pow(lambda,6)) + 280*pow(t,4)*(2*pow(r0,4)*pow(lambda,4) - 9*pow(r0,2)*pow(lambda,6) + 3*pow(lambda,8)) - 
           28*pow(t,3)*(12*pow(r0,5)*pow(lambda,5) - 100*pow(r0,3)*pow(lambda,7) + 105*r0*pow(lambda,9))) + 
        pow(E,pow(t - r0*lambda,2)/pow(lambda,4))*(16*pow(t,8) + 112*r0*pow(t,7)*lambda + 112*pow(t,6)*(3*pow(r0,2)*pow(lambda,2) - 2*pow(lambda,4)) + 56*pow(t,5)*(10*pow(r0,3)*pow(lambda,3) - 21*r0*pow(lambda,5)) + 
           2*r0*t*pow(lambda,7)*(8*pow(r0,6) - 252*pow(r0,4)*pow(lambda,2) + 1050*pow(r0,2)*pow(lambda,4) - 735*pow(lambda,6)) + 
           28*pow(t,2)*pow(lambda,6)*(4*pow(r0,6) - 60*pow(r0,4)*pow(lambda,2) + 135*pow(r0,2)*pow(lambda,4) - 30*pow(lambda,6)) + 7*pow(lambda,10)*(-8*pow(r0,6) + 60*pow(r0,4)*pow(lambda,2) - 90*pow(r0,2)*pow(lambda,4) + 15*pow(lambda,6)) + 
           280*pow(t,4)*(2*pow(r0,4)*pow(lambda,4) - 9*pow(r0,2)*pow(lambda,6) + 3*pow(lambda,8)) + 28*pow(t,3)*(12*pow(r0,5)*pow(lambda,5) - 100*pow(r0,3)*pow(lambda,7) + 105*r0*pow(lambda,9)))))/
    (105.*pow(E,(2*(pow(t,2) + pow(r0,2)*pow(lambda,2)))/pow(lambda,4))*pow(lambda,24)) + (A*pow(rr,6)*(-1 + 11*lambda*(1 + 5*lambda*(-1 + 3*lambda)))*
      (pow(E,pow(t - r0*lambda,2)/pow(lambda,4))*(64*pow(t,12) + 704*r0*pow(t,11)*lambda + 704*pow(t,10)*(5*pow(r0,2)*pow(lambda,2) - 3*pow(lambda,4)) + 1760*pow(t,9)*(6*pow(r0,3)*pow(lambda,3) - 11*r0*pow(lambda,5)) + 
           528*r0*pow(t,5)*pow(lambda,7)*(40*pow(r0,6) - 588*pow(r0,4)*pow(lambda,2) + 1890*pow(r0,2)*pow(lambda,4) - 1155*pow(lambda,6)) + 
           7392*pow(t,6)*pow(lambda,6)*(4*pow(r0,6) - 40*pow(r0,4)*pow(lambda,2) + 75*pow(r0,2)*pow(lambda,4) - 15*pow(lambda,6)) + 2640*pow(t,8)*(8*pow(r0,4)*pow(lambda,4) - 30*pow(r0,2)*pow(lambda,6) + 9*pow(lambda,8)) + 
           660*pow(t,4)*pow(lambda,8)*(16*pow(r0,8) - 336*pow(r0,6)*pow(lambda,2) + 1680*pow(r0,4)*pow(lambda,4) - 2100*pow(r0,2)*pow(lambda,6) + 315*pow(lambda,8)) + 
           220*r0*pow(t,3)*pow(lambda,9)*(16*pow(r0,8) - 480*pow(r0,6)*pow(lambda,2) + 3528*pow(r0,4)*pow(lambda,4) - 7560*pow(r0,2)*pow(lambda,6) + 3465*pow(lambda,8)) + 
           1056*pow(t,7)*(28*pow(r0,5)*pow(lambda,5) - 180*pow(r0,3)*pow(lambda,7) + 165*r0*pow(lambda,9)) + 
           2*r0*t*pow(lambda,11)*(32*pow(r0,10) - 2640*pow(r0,8)*pow(lambda,2) + 39600*pow(r0,6)*pow(lambda,4) - 194040*pow(r0,4)*pow(lambda,6) + 311850*pow(r0,2)*pow(lambda,8) - 114345*pow(lambda,10)) + 
           44*pow(t,2)*pow(lambda,10)*(16*pow(r0,10) - 720*pow(r0,8)*pow(lambda,2) + 7560*pow(r0,6)*pow(lambda,4) - 25200*pow(r0,4)*pow(lambda,6) + 23625*pow(r0,2)*pow(lambda,8) - 2835*pow(lambda,10)) + 
           11*pow(lambda,14)*(-32*pow(r0,10) + 720*pow(r0,8)*pow(lambda,2) - 5040*pow(r0,6)*pow(lambda,4) + 12600*pow(r0,4)*pow(lambda,6) - 9450*pow(r0,2)*pow(lambda,8) + 945*pow(lambda,10))) + 
        pow(E,pow(t + r0*lambda,2)/pow(lambda,4))*(64*pow(t,12) - 704*r0*pow(t,11)*lambda + 704*pow(t,10)*(5*pow(r0,2)*pow(lambda,2) - 3*pow(lambda,4)) - 1760*pow(t,9)*(6*pow(r0,3)*pow(lambda,3) - 11*r0*pow(lambda,5)) - 
           528*r0*pow(t,5)*pow(lambda,7)*(40*pow(r0,6) - 588*pow(r0,4)*pow(lambda,2) + 1890*pow(r0,2)*pow(lambda,4) - 1155*pow(lambda,6)) + 
           7392*pow(t,6)*pow(lambda,6)*(4*pow(r0,6) - 40*pow(r0,4)*pow(lambda,2) + 75*pow(r0,2)*pow(lambda,4) - 15*pow(lambda,6)) + 2640*pow(t,8)*(8*pow(r0,4)*pow(lambda,4) - 30*pow(r0,2)*pow(lambda,6) + 9*pow(lambda,8)) + 
           660*pow(t,4)*pow(lambda,8)*(16*pow(r0,8) - 336*pow(r0,6)*pow(lambda,2) + 1680*pow(r0,4)*pow(lambda,4) - 2100*pow(r0,2)*pow(lambda,6) + 315*pow(lambda,8)) - 
           220*r0*pow(t,3)*pow(lambda,9)*(16*pow(r0,8) - 480*pow(r0,6)*pow(lambda,2) + 3528*pow(r0,4)*pow(lambda,4) - 7560*pow(r0,2)*pow(lambda,6) + 3465*pow(lambda,8)) - 
           1056*pow(t,7)*(28*pow(r0,5)*pow(lambda,5) - 180*pow(r0,3)*pow(lambda,7) + 165*r0*pow(lambda,9)) + 
           44*pow(t,2)*pow(lambda,10)*(16*pow(r0,10) - 720*pow(r0,8)*pow(lambda,2) + 7560*pow(r0,6)*pow(lambda,4) - 25200*pow(r0,4)*pow(lambda,6) + 23625*pow(r0,2)*pow(lambda,8) - 2835*pow(lambda,10)) + 
           11*pow(lambda,14)*(-32*pow(r0,10) + 720*pow(r0,8)*pow(lambda,2) - 5040*pow(r0,6)*pow(lambda,4) + 12600*pow(r0,4)*pow(lambda,6) - 9450*pow(r0,2)*pow(lambda,8) + 945*pow(lambda,10)) + 
           2*r0*t*pow(lambda,11)*(-32*pow(r0,10) + 2640*pow(r0,8)*pow(lambda,2) - 39600*pow(r0,6)*pow(lambda,4) + 194040*pow(r0,4)*pow(lambda,6) - 311850*pow(r0,2)*pow(lambda,8) + 114345*pow(lambda,10)))))/
    (51975.*pow(E,(2*(pow(t,2) + pow(r0,2)*pow(lambda,2)))/pow(lambda,4))*pow(lambda,40));

C2 = (6*A*(-1 + lambda)*(-2*pow(t,2) - 2*r0*t*lambda + pow(lambda,4) + pow(E,(4*r0*t)/pow(lambda,3))*(-2*pow(t,2) + 2*r0*t*lambda + pow(lambda,4))))/(pow(E,pow(t + r0*lambda,2)/pow(lambda,4))*pow(rr,4)) - 
   (2*A*(-1 + lambda)*(1 - 2*lambda + 4*pow(lambda,2))*(pow(E,pow(t - r0*lambda,2)/pow(lambda,4))*(4*pow(t,4) + 12*r0*pow(t,3)*lambda - 6*pow(r0,2)*pow(lambda,6) + 3*pow(lambda,8) + 12*pow(t,2)*(r0 - lambda)*pow(lambda,2)*(r0 + lambda) + 
           2*t*(2*pow(r0,3)*pow(lambda,3) - 9*r0*pow(lambda,5))) + pow(E,pow(t + r0*lambda,2)/pow(lambda,4))*
         (4*pow(t,4) - 12*r0*pow(t,3)*lambda - 6*pow(r0,2)*pow(lambda,6) + 3*pow(lambda,8) + 12*pow(t,2)*(r0 - lambda)*pow(lambda,2)*(r0 + lambda) + t*(-4*pow(r0,3)*pow(lambda,3) + 18*r0*pow(lambda,5)))))/
    (pow(E,(2*(pow(t,2) + pow(r0,2)*pow(lambda,2)))/pow(lambda,4))*pow(rr,2)*pow(lambda,8)) - (A*(1 + 5*lambda*(-1 + 4*lambda*(1 + 2*(-1 + lambda)*lambda)))*
      (pow(E,pow(t - r0*lambda,2)/pow(lambda,4))*(-8*pow(t,6) - 40*r0*pow(t,5)*lambda + 40*r0*pow(t,3)*pow(lambda,3)*(-2*pow(r0,2) + 5*pow(lambda,2)) + 5*pow(lambda,8)*(4*pow(r0,4) - 12*pow(r0,2)*pow(lambda,2) + 3*pow(lambda,4)) + 
           pow(t,4)*(-80*pow(r0,2)*pow(lambda,2) + 60*pow(lambda,4)) - 10*pow(t,2)*(4*pow(r0,4)*pow(lambda,4) - 24*pow(r0,2)*pow(lambda,6) + 9*pow(lambda,8)) - 2*t*(4*pow(r0,5)*pow(lambda,5) - 60*pow(r0,3)*pow(lambda,7) + 75*r0*pow(lambda,9))) + 
        pow(E,pow(t + r0*lambda,2)/pow(lambda,4))*(-8*pow(t,6) + 40*r0*pow(t,5)*lambda + 5*pow(lambda,8)*(4*pow(r0,4) - 12*pow(r0,2)*pow(lambda,2) + 3*pow(lambda,4)) + pow(t,4)*(-80*pow(r0,2)*pow(lambda,2) + 60*pow(lambda,4)) + 
           40*pow(t,3)*(2*pow(r0,3)*pow(lambda,3) - 5*r0*pow(lambda,5)) - 10*pow(t,2)*(4*pow(r0,4)*pow(lambda,4) - 24*pow(r0,2)*pow(lambda,6) + 9*pow(lambda,8)) + 2*t*(4*pow(r0,5)*pow(lambda,5) - 60*pow(r0,3)*pow(lambda,7) + 75*r0*pow(lambda,9)))
        ))/(5.*pow(E,(2*(pow(t,2) + pow(r0,2)*pow(lambda,2)))/pow(lambda,4))*pow(lambda,16)) + (A*pow(rr,2)*
      (pow(E,pow(t + r0*lambda,2)/pow(lambda,4))*(16*pow(t,8) - 112*r0*pow(t,7)*lambda + 56*r0*pow(t,5)*pow(lambda,3)*(-10*pow(r0,2) + 21*pow(lambda,2)) + 112*pow(t,6)*(3*pow(r0,2)*pow(lambda,2) - 2*pow(lambda,4)) + 
           28*pow(t,2)*pow(lambda,6)*(4*pow(r0,6) - 60*pow(r0,4)*pow(lambda,2) + 135*pow(r0,2)*pow(lambda,4) - 30*pow(lambda,6)) + 7*pow(lambda,10)*(-8*pow(r0,6) + 60*pow(r0,4)*pow(lambda,2) - 90*pow(r0,2)*pow(lambda,4) + 15*pow(lambda,6)) + 
           2*r0*t*pow(lambda,7)*(-8*pow(r0,6) + 252*pow(r0,4)*pow(lambda,2) - 1050*pow(r0,2)*pow(lambda,4) + 735*pow(lambda,6)) + 280*pow(t,4)*(2*pow(r0,4)*pow(lambda,4) - 9*pow(r0,2)*pow(lambda,6) + 3*pow(lambda,8)) - 
           28*pow(t,3)*(12*pow(r0,5)*pow(lambda,5) - 100*pow(r0,3)*pow(lambda,7) + 105*r0*pow(lambda,9))) + 
        pow(E,pow(t - r0*lambda,2)/pow(lambda,4))*(16*pow(t,8) + 112*r0*pow(t,7)*lambda + 112*pow(t,6)*(3*pow(r0,2)*pow(lambda,2) - 2*pow(lambda,4)) + 56*pow(t,5)*(10*pow(r0,3)*pow(lambda,3) - 21*r0*pow(lambda,5)) + 
           2*r0*t*pow(lambda,7)*(8*pow(r0,6) - 252*pow(r0,4)*pow(lambda,2) + 1050*pow(r0,2)*pow(lambda,4) - 735*pow(lambda,6)) + 
           28*pow(t,2)*pow(lambda,6)*(4*pow(r0,6) - 60*pow(r0,4)*pow(lambda,2) + 135*pow(r0,2)*pow(lambda,4) - 30*pow(lambda,6)) + 7*pow(lambda,10)*(-8*pow(r0,6) + 60*pow(r0,4)*pow(lambda,2) - 90*pow(r0,2)*pow(lambda,4) + 15*pow(lambda,6)) + 
           280*pow(t,4)*(2*pow(r0,4)*pow(lambda,4) - 9*pow(r0,2)*pow(lambda,6) + 3*pow(lambda,8)) + 28*pow(t,3)*(12*pow(r0,5)*pow(lambda,5) - 100*pow(r0,3)*pow(lambda,7) + 105*r0*pow(lambda,9))))*(1 + 7*lambda*(-1 + 2*lambda*(3 + 10*lambda*(-1 + 2*lambda)))))/
    (105.*pow(E,(2*(pow(t,2) + pow(r0,2)*pow(lambda,2)))/pow(lambda,4))*pow(lambda,24)) - (A*pow(rr,4)*
      (pow(E,pow(t - r0*lambda,2)/pow(lambda,4))*(-32*pow(t,10) - 288*r0*pow(t,9)*lambda + 144*pow(t,8)*pow(lambda,2)*(-8*pow(r0,2) + 5*pow(lambda,2)) + 192*r0*pow(t,7)*pow(lambda,3)*(-14*pow(r0,2) + 27*pow(lambda,2)) - 
           144*r0*pow(t,3)*pow(lambda,7)*(8*pow(r0,6) - 140*pow(r0,4)*pow(lambda,2) + 490*pow(r0,2)*pow(lambda,4) - 315*pow(lambda,6)) + 
           168*pow(t,4)*pow(lambda,6)*(-16*pow(r0,6) + 180*pow(r0,4)*pow(lambda,2) - 360*pow(r0,2)*pow(lambda,4) + 75*pow(lambda,6)) - 1008*pow(t,6)*(4*pow(r0,4)*pow(lambda,4) - 16*pow(r0,2)*pow(lambda,6) + 5*pow(lambda,8)) + 
           9*pow(lambda,12)*(16*pow(r0,8) - 224*pow(r0,6)*pow(lambda,2) + 840*pow(r0,4)*pow(lambda,4) - 840*pow(r0,2)*pow(lambda,6) + 105*pow(lambda,8)) - 
           18*pow(t,2)*pow(lambda,8)*(16*pow(r0,8) - 448*pow(r0,6)*pow(lambda,2) + 2520*pow(r0,4)*pow(lambda,4) - 3360*pow(r0,2)*pow(lambda,6) + 525*pow(lambda,8)) - 
           2*r0*t*pow(lambda,9)*(16*pow(r0,8) - 864*pow(r0,6)*pow(lambda,2) + 7560*pow(r0,4)*pow(lambda,4) - 17640*pow(r0,2)*pow(lambda,6) + 8505*pow(lambda,8)) - 
           1008*pow(t,5)*(4*pow(r0,5)*pow(lambda,5) - 28*pow(r0,3)*pow(lambda,7) + 27*r0*pow(lambda,9))) + 
        pow(E,pow(t + r0*lambda,2)/pow(lambda,4))*(-32*pow(t,10) + 288*r0*pow(t,9)*lambda + 144*pow(t,8)*pow(lambda,2)*(-8*pow(r0,2) + 5*pow(lambda,2)) + 192*pow(t,7)*(14*pow(r0,3)*pow(lambda,3) - 27*r0*pow(lambda,5)) + 
           144*r0*pow(t,3)*pow(lambda,7)*(8*pow(r0,6) - 140*pow(r0,4)*pow(lambda,2) + 490*pow(r0,2)*pow(lambda,4) - 315*pow(lambda,6)) + 
           168*pow(t,4)*pow(lambda,6)*(-16*pow(r0,6) + 180*pow(r0,4)*pow(lambda,2) - 360*pow(r0,2)*pow(lambda,4) + 75*pow(lambda,6)) - 1008*pow(t,6)*(4*pow(r0,4)*pow(lambda,4) - 16*pow(r0,2)*pow(lambda,6) + 5*pow(lambda,8)) + 
           9*pow(lambda,12)*(16*pow(r0,8) - 224*pow(r0,6)*pow(lambda,2) + 840*pow(r0,4)*pow(lambda,4) - 840*pow(r0,2)*pow(lambda,6) + 105*pow(lambda,8)) - 
           18*pow(t,2)*pow(lambda,8)*(16*pow(r0,8) - 448*pow(r0,6)*pow(lambda,2) + 2520*pow(r0,4)*pow(lambda,4) - 3360*pow(r0,2)*pow(lambda,6) + 525*pow(lambda,8)) + 
           2*r0*t*pow(lambda,9)*(16*pow(r0,8) - 864*pow(r0,6)*pow(lambda,2) + 7560*pow(r0,4)*pow(lambda,4) - 17640*pow(r0,2)*pow(lambda,6) + 8505*pow(lambda,8)) + 
           1008*pow(t,5)*(4*pow(r0,5)*pow(lambda,5) - 28*pow(r0,3)*pow(lambda,7) + 27*r0*pow(lambda,9))))*(1 + 3*lambda*(-3 + 8*lambda*(3 + 14*lambda*(-1 + 3*lambda)))))/(3780.*pow(E,(2*(pow(t,2) + pow(r0,2)*pow(lambda,2)))/pow(lambda,4))*pow(lambda,32)) + 
   (A*pow(rr,6)*(pow(E,pow(t - r0*lambda,2)/pow(lambda,4))*(64*pow(t,12) + 704*r0*pow(t,11)*lambda + 704*pow(t,10)*(5*pow(r0,2)*pow(lambda,2) - 3*pow(lambda,4)) + 1760*pow(t,9)*(6*pow(r0,3)*pow(lambda,3) - 11*r0*pow(lambda,5)) + 
           528*r0*pow(t,5)*pow(lambda,7)*(40*pow(r0,6) - 588*pow(r0,4)*pow(lambda,2) + 1890*pow(r0,2)*pow(lambda,4) - 1155*pow(lambda,6)) + 
           7392*pow(t,6)*pow(lambda,6)*(4*pow(r0,6) - 40*pow(r0,4)*pow(lambda,2) + 75*pow(r0,2)*pow(lambda,4) - 15*pow(lambda,6)) + 2640*pow(t,8)*(8*pow(r0,4)*pow(lambda,4) - 30*pow(r0,2)*pow(lambda,6) + 9*pow(lambda,8)) + 
           660*pow(t,4)*pow(lambda,8)*(16*pow(r0,8) - 336*pow(r0,6)*pow(lambda,2) + 1680*pow(r0,4)*pow(lambda,4) - 2100*pow(r0,2)*pow(lambda,6) + 315*pow(lambda,8)) + 
           220*r0*pow(t,3)*pow(lambda,9)*(16*pow(r0,8) - 480*pow(r0,6)*pow(lambda,2) + 3528*pow(r0,4)*pow(lambda,4) - 7560*pow(r0,2)*pow(lambda,6) + 3465*pow(lambda,8)) + 
           1056*pow(t,7)*(28*pow(r0,5)*pow(lambda,5) - 180*pow(r0,3)*pow(lambda,7) + 165*r0*pow(lambda,9)) + 
           2*r0*t*pow(lambda,11)*(32*pow(r0,10) - 2640*pow(r0,8)*pow(lambda,2) + 39600*pow(r0,6)*pow(lambda,4) - 194040*pow(r0,4)*pow(lambda,6) + 311850*pow(r0,2)*pow(lambda,8) - 114345*pow(lambda,10)) + 
           44*pow(t,2)*pow(lambda,10)*(16*pow(r0,10) - 720*pow(r0,8)*pow(lambda,2) + 7560*pow(r0,6)*pow(lambda,4) - 25200*pow(r0,4)*pow(lambda,6) + 23625*pow(r0,2)*pow(lambda,8) - 2835*pow(lambda,10)) + 
           11*pow(lambda,14)*(-32*pow(r0,10) + 720*pow(r0,8)*pow(lambda,2) - 5040*pow(r0,6)*pow(lambda,4) + 12600*pow(r0,4)*pow(lambda,6) - 9450*pow(r0,2)*pow(lambda,8) + 945*pow(lambda,10))) + 
        pow(E,pow(t + r0*lambda,2)/pow(lambda,4))*(64*pow(t,12) - 704*r0*pow(t,11)*lambda + 704*pow(t,10)*(5*pow(r0,2)*pow(lambda,2) - 3*pow(lambda,4)) - 1760*pow(t,9)*(6*pow(r0,3)*pow(lambda,3) - 11*r0*pow(lambda,5)) - 
           528*r0*pow(t,5)*pow(lambda,7)*(40*pow(r0,6) - 588*pow(r0,4)*pow(lambda,2) + 1890*pow(r0,2)*pow(lambda,4) - 1155*pow(lambda,6)) + 
           7392*pow(t,6)*pow(lambda,6)*(4*pow(r0,6) - 40*pow(r0,4)*pow(lambda,2) + 75*pow(r0,2)*pow(lambda,4) - 15*pow(lambda,6)) + 2640*pow(t,8)*(8*pow(r0,4)*pow(lambda,4) - 30*pow(r0,2)*pow(lambda,6) + 9*pow(lambda,8)) + 
           660*pow(t,4)*pow(lambda,8)*(16*pow(r0,8) - 336*pow(r0,6)*pow(lambda,2) + 1680*pow(r0,4)*pow(lambda,4) - 2100*pow(r0,2)*pow(lambda,6) + 315*pow(lambda,8)) - 
           220*r0*pow(t,3)*pow(lambda,9)*(16*pow(r0,8) - 480*pow(r0,6)*pow(lambda,2) + 3528*pow(r0,4)*pow(lambda,4) - 7560*pow(r0,2)*pow(lambda,6) + 3465*pow(lambda,8)) - 
           1056*pow(t,7)*(28*pow(r0,5)*pow(lambda,5) - 180*pow(r0,3)*pow(lambda,7) + 165*r0*pow(lambda,9)) + 
           44*pow(t,2)*pow(lambda,10)*(16*pow(r0,10) - 720*pow(r0,8)*pow(lambda,2) + 7560*pow(r0,6)*pow(lambda,4) - 25200*pow(r0,4)*pow(lambda,6) + 23625*pow(r0,2)*pow(lambda,8) - 2835*pow(lambda,10)) + 
           11*pow(lambda,14)*(-32*pow(r0,10) + 720*pow(r0,8)*pow(lambda,2) - 5040*pow(r0,6)*pow(lambda,4) + 12600*pow(r0,4)*pow(lambda,6) - 9450*pow(r0,2)*pow(lambda,8) + 945*pow(lambda,10)) + 
           2*r0*t*pow(lambda,11)*(-32*pow(r0,10) + 2640*pow(r0,8)*pow(lambda,2) - 39600*pow(r0,6)*pow(lambda,4) + 194040*pow(r0,4)*pow(lambda,6) - 311850*pow(r0,2)*pow(lambda,8) + 114345*pow(lambda,10))))*
      (1 + 11*lambda*(-1 + 10*lambda*(1 + 6*lambda*(-1 + 4*lambda)))))/(207900.*pow(E,(2*(pow(t,2) + pow(r0,2)*pow(lambda,2)))/pow(lambda,4))*pow(lambda,40));
//This is here to safeguard when t gets high at the origin we have exp (-400) These terms just go to 0
// After computing A2, B2, C2
if (std::isnan(A2) || std::isinf(A2)) A2 = 0.0;
if (std::isnan(B2) || std::isinf(B2)) B2 = 0.0;
if (std::isnan(C2) || std::isinf(C2)) C2 = 0.0;
    }
//else we continue with the numerical formulation:
else{

A2 = (12*A*((3*(t + rr*(-1 + lambda))*pow(lambda,6) - 4*pow(rr,2)*(rr - t)*pow(-rr + t + r0*lambda,2) - 2*rr*pow(lambda,3)*(3*pow(rr - t,2) - 3*(rr + r0)*(rr - t)*lambda + 2*rr*r0*pow(lambda,2)))/pow(E,pow(-rr + t + r0*lambda,2)/pow(lambda,4)) - 
       (3*pow(lambda,6)*(rr + t - rr*lambda) + 4*pow(rr,2)*(rr + t)*pow(rr + t - r0*lambda,2) + 2*rr*pow(lambda,3)*(3*pow(rr + t,2) - 3*(rr + r0)*(rr + t)*lambda + 2*rr*r0*pow(lambda,2)))/pow(E,pow(rr + t - r0*lambda,2)/pow(lambda,4)) + 
       (3*(t + rr*(-1 + lambda))*pow(lambda,6) - 4*pow(rr,2)*(rr - t)*pow(rr - t + r0*lambda,2) - 2*rr*pow(lambda,3)*(rr*lambda*(-3*rr + 3*t - 2*r0*lambda) + 3*(rr - t)*(rr - t + r0*lambda)))/pow(E,pow(rr - t + r0*lambda,2)/pow(lambda,4)) - 
       (3*pow(lambda,6)*(rr + t - rr*lambda) + 4*pow(rr,2)*(rr + t)*pow(rr + t + r0*lambda,2) + 2*rr*pow(lambda,3)*(3*(rr + t)*(rr + t + r0*lambda) - rr*lambda*(3*(rr + t) + 2*r0*lambda)))/pow(E,pow(rr + t + r0*lambda,2)/pow(lambda,4))))/(pow(rr,5)*pow(lambda,2));

B2 = (4*A*(3*((pow(E,-pow(rr - t + r0*lambda,2)/pow(lambda,4)) + pow(E,-pow(-rr + t + r0*lambda,2)/pow(lambda,4)))*(rr - t) + (pow(E,-pow(rr + t - r0*lambda,2)/pow(lambda,4)) + pow(E,-pow(rr + t + r0*lambda,2)/pow(lambda,4)))*(rr + t))*pow(lambda,9) - 
       3*rr*pow(lambda,6)*((pow(lambda,4) - 2*(rr - t)*(rr - t - r0*lambda))/pow(E,pow(-rr + t + r0*lambda,2)/pow(lambda,4)) + (pow(lambda,4) - 2*(rr + t)*(rr + t - r0*lambda))/pow(E,pow(rr + t - r0*lambda,2)/pow(lambda,4)) + 
          (pow(lambda,4) - 2*(rr - t)*(rr - t + r0*lambda))/pow(E,pow(rr - t + r0*lambda,2)/pow(lambda,4)) + (pow(lambda,4) - 2*(rr + t)*(rr + t + r0*lambda))/pow(E,pow(rr + t + r0*lambda,2)/pow(lambda,4))) + 
       pow(rr,3)*((3*pow(lambda,8) + 4*(rr + t)*pow(rr + t - r0*lambda,3) - 6*pow(lambda,4)*(rr + t - r0*lambda)*(2*(rr + t) - r0*lambda))/pow(E,pow(rr + t - r0*lambda,2)/pow(lambda,4)) + 
          (3*pow(lambda,8) - 6*pow(lambda,4)*(2*rr - 2*t + r0*lambda)*(rr - t + r0*lambda) + 4*(rr - t)*pow(rr - t + r0*lambda,3))/pow(E,pow(rr - t + r0*lambda,2)/pow(lambda,4)) + 
          (3*pow(lambda,8) + 4*(rr - t)*pow(rr - t - r0*lambda,3) - 6*pow(lambda,4)*(-rr + t + r0*lambda)*(-2*rr + 2*t + r0*lambda))/pow(E,pow(-rr + t + r0*lambda,2)/pow(lambda,4)) + 
          (3*pow(lambda,8) + 4*(rr + t)*pow(rr + t + r0*lambda,3) - 6*pow(lambda,4)*(rr + t + r0*lambda)*(2*(rr + t) + r0*lambda))/pow(E,pow(rr + t + r0*lambda,2)/pow(lambda,4))) - 
       3*pow(rr,2)*pow(lambda,3)*((pow(lambda,4)*(3*(rr + t) - 2*r0*lambda) - 2*(rr + t)*pow(rr + t - r0*lambda,2))/pow(E,pow(rr + t - r0*lambda,2)/pow(lambda,4)) + 
          (pow(lambda,4)*(3*rr - 3*t - 2*r0*lambda) - 2*(rr - t)*pow(-rr + t + r0*lambda,2))/pow(E,pow(-rr + t + r0*lambda,2)/pow(lambda,4)) + 
          (-2*(rr - t)*pow(rr - t + r0*lambda,2) + pow(lambda,4)*(3*rr - 3*t + 2*r0*lambda))/pow(E,pow(rr - t + r0*lambda,2)/pow(lambda,4)) + (-2*(rr + t)*pow(rr + t + r0*lambda,2) + pow(lambda,4)*(3*(rr + t) + 2*r0*lambda))/pow(E,pow(rr + t + r0*lambda,2)/pow(lambda,4)))))
    /(pow(rr,5)*pow(lambda,5));

C2 =-((A*((3*pow(lambda,12)*(rr + t - rr*lambda) + 16*pow(rr,4)*(rr + t)*pow(rr + t - r0*lambda,4) + 6*rr*pow(lambda,9)*(rr + t - 2*rr*lambda)*(rr + t - (rr + r0)*lambda) + 
            16*pow(rr,3)*pow(lambda,3)*pow(rr + t - r0*lambda,2)*(pow(rr + t,2) - (5*rr + r0)*(rr + t)*lambda + 2*rr*r0*pow(lambda,2)) + 
            12*pow(rr,2)*pow(lambda,6)*(pow(rr + t,3) - 2*(2*rr + r0)*pow(rr + t,2)*lambda + (rr + r0)*(5*rr + r0)*(rr + t)*pow(lambda,2) - 2*rr*r0*(2*rr + r0)*pow(lambda,3)))/pow(E,pow(rr + t - r0*lambda,2)/pow(lambda,4)) + 
         (-3*(t + rr*(-1 + lambda))*pow(lambda,12) + 16*pow(rr,4)*(rr - t)*pow(-rr + t + r0*lambda,4) + 16*pow(rr,3)*pow(lambda,3)*pow(-rr + t + r0*lambda,2)*(pow(rr - t,2) - (5*rr + r0)*(rr - t)*lambda + 2*rr*r0*pow(lambda,2)) + 
            12*pow(rr,2)*pow(lambda,6)*(pow(rr - t,3) - 2*(2*rr + r0)*pow(rr - t,2)*lambda + (rr + r0)*(5*rr + r0)*(rr - t)*pow(lambda,2) - 2*rr*r0*(2*rr + r0)*pow(lambda,3)) + 6*rr*pow(lambda,9)*(t + rr*(-1 + lambda) + r0*lambda)*(t + rr*(-1 + 2*lambda)))/
          pow(E,pow(-rr + t + r0*lambda,2)/pow(lambda,4)) + (-3*(t + rr*(-1 + lambda))*pow(lambda,12) + 16*pow(rr,4)*(rr - t)*pow(rr - t + r0*lambda,4) + 6*rr*pow(lambda,9)*(t + rr*(-1 + lambda) - r0*lambda)*(t + rr*(-1 + 2*lambda)) + 
            16*pow(rr,3)*pow(lambda,3)*pow(rr - t + r0*lambda,2)*(rr*lambda*(-5*rr + 5*t - 2*r0*lambda) + (rr - t)*(rr - t + r0*lambda)) - 
            12*pow(rr,2)*pow(lambda,6)*(pow(rr,2)*pow(lambda,2)*(-5*rr + 5*t - 4*r0*lambda) + 2*rr*lambda*(2*rr - 2*t + r0*lambda)*(rr - t + r0*lambda) - (rr - t)*pow(rr - t + r0*lambda,2)))/pow(E,pow(rr - t + r0*lambda,2)/pow(lambda,4)) + 
         (3*pow(lambda,12)*(rr + t - rr*lambda) + 16*pow(rr,4)*(rr + t)*pow(rr + t + r0*lambda,4) + 6*rr*pow(lambda,9)*(rr + t - 2*rr*lambda)*(rr + t - rr*lambda + r0*lambda) + 
            16*pow(rr,3)*pow(lambda,3)*pow(rr + t + r0*lambda,2)*((rr + t)*(rr + t + r0*lambda) - rr*lambda*(5*(rr + t) + 2*r0*lambda)) + 
            12*pow(rr,2)*pow(lambda,6)*((rr + t)*pow(rr + t + r0*lambda,2) - 2*rr*lambda*(rr + t + r0*lambda)*(2*(rr + t) + r0*lambda) + pow(rr,2)*pow(lambda,2)*(5*(rr + t) + 4*r0*lambda)))/pow(E,pow(rr + t + r0*lambda,2)/pow(lambda,4))))/(pow(rr,5)*pow(lambda,8)));
}

// double Y20= 2- 3* sin_theta*sin_theta;
// double Y20theta =-6*cos_theta*sin_theta;
// double Y20_thetatheta = 3* sin_theta*sin_theta;

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

// Now for the angular functions for l=2
///////////////////////////////////////////////////////////////////////////
// Starting with Y_20
//
Y20        = (1.0 / 4.0) * sqrt(5.0 / pi) * (3.0 * pow(cos_theta, 2) - 1.0);
Y21        = (1.0 / 2.0) * sqrt(15.0 / pi) * sin_theta * cos_theta * cos_phi;
Y2n1       = (1.0 / 2.0) * sqrt(15.0 / pi) * sin_theta * cos_theta * sin_phi;
Y22        = (1.0 / 4.0) * sqrt(15.0 / pi) * pow(sin_theta, 2) * cos_2phi;
Y2n2       = (1.0 / 4.0) * sqrt(15.0 / pi) * pow(sin_theta, 2) * sin_2phi;

// Y_2m,theta
//
Y20_theta  = -(3.0 / 2.0) * sqrt(5.0 / pi) * sin_theta * cos_theta;
Y21_theta  = (1.0 / 2.0) * sqrt(15.0 / pi) * cos_2theta * cos_phi;
Y2n1_theta = (1.0 / 2.0) * sqrt(15.0 / pi) * cos_2theta * sin_phi;
Y22_theta  = (1.0 / 2.0) * sqrt(15.0 / pi) * sin_theta * cos_theta * cos_2phi;
Y2n2_theta = (1.0 / 2.0) * sqrt(15.0 / pi) * sin_theta * cos_theta * sin_2phi;

// Y_2m,phi
//
Y20_phi    = 0.0;
Y21_phi    = -(1.0 / 2.0) * sqrt(15.0 / pi) * cos_theta * sin_phi;
Y2n1_phi   = (1.0 / 2.0) * sqrt(15.0 / pi) * cos_theta * cos_phi;
Y22_phi    = -(1.0 / 2.0) * sqrt(15.0 / pi) * sin_theta * sin_2phi;
Y2n2_phi   = (1.0 / 2.0) * sqrt(15.0 / pi) * sin_theta * cos_2phi;

// Y_2m,theta,theta + .5 l ( l + 1 ) Y_2m
//
Y20_thth   = (3.0 / 4.0) * sqrt(5.0 / pi) * pow(sin_theta, 2);
Y21_thth   = -(1.0 / 2.0) * sqrt(15.0 / pi) * sin_theta * cos_theta * cos_phi;
Y2n1_thth  = -(1.0 / 2.0) * sqrt(15.0 / pi) * sin_theta * cos_theta * sin_phi;
Y22_thth   = (1.0 / 8.0) * sqrt(15.0 / pi) * (3.0 + cos_2theta) * cos_2phi;
Y2n2_thth  = (1.0 / 8.0) * sqrt(15.0 / pi) * (3.0 + cos_2theta) * sin_2phi;

// Y_2m,theta,phi - cot_theta Y-2m,phi
//
Y20_thphi  = 0.0;
Y21_thphi  = (1.0 / 2.0) * sqrt(15.0 / pi) * sin_theta * sin_phi;
Y2n1_thphi = -(1.0 / 2.0) * sqrt(15.0 / pi) * sin_theta * cos_phi;
Y22_thphi  = -(1.0 / 2.0) * sqrt(15.0 / pi) * cos_theta * sin_2phi;
Y2n2_thphi = (1.0 / 2.0) * sqrt(15.0 / pi) * cos_theta * cos_2phi;

        Ylm       = Y20;
        Ylm_theta = Y20_theta;
        Ylm_thth  = Y20_thth;
        Ylm_thphi = Y20_thphi;
        Ylm_phi   = Y20_phi;

//converting into cartesian coordinates.....

//defining some things for the conversion://
// double g_bare_rr   = g_rr-1.0;
// double g_bare_thth = g_thth -rr*rr;
// double g_bare_phph = g_phph  - pow(rr,2)*pow(sin_theta,2);

// double h_bare_rr   = A2* Y20;
// double h_bare_rth  = (2*B2*Y20theta);
// double h_bare_rph  = 0.0;
// double h_bare_thth = A2*Y20/2.0 + C2*Y20_thetatheta;
// double h_bare_thph = 0.0;
// double h_bare_phph = A2*Y20/2.0 +  C2*Y20_thetatheta;

// double g_bare_rr   = 1.0 + h_bare_rr;
// double g_bare_rth  = h_bare_rth;
// double g_bare_rph  = h_bare_rph;
// double g_bare_thth = 1.0 - h_bare_thth;
// double g_bare_thph = h_bare_thph;
// double g_bare_phph = 1.0 - h_bare_phph;

// double g_rr        = g_bare_rr;
// double g_rth       = rr * g_bare_rth;
// double g_rph       = rr * sin_theta * g_bare_rph;
// double g_thth      = rr * rr * g_bare_thth;
// double g_thph      = rr * rr * sin_theta * g_bare_thph;
// double g_phph      = rr * rr * sin_theta * sin_theta * g_bare_phph;
////////////////////////////////////////////////////////////////////////////////
// Now put it all together to get the perturbed metric

A_lm= A2;
B_lm= B2;
C_lm= C2;

h_bare_rr   = A_lm * Ylm;
h_bare_rth  = B_lm * Ylm_theta - K_lm * Ylm_phi;
h_bare_rph  = B_lm * Ylm_phi + K_lm * Ylm_theta;
h_bare_thth = -0.5 * A_lm * Ylm + C_lm * Ylm_thth - L_lm * Ylm_thphi;
h_bare_thph = C_lm * Ylm_thphi + L_lm * Ylm_thth;
h_bare_phph = -(0.5 * A_lm * Ylm + C_lm * Ylm_thth - L_lm * Ylm_thphi);

g_bare_rr   = 1.0 + h_bare_rr;
g_bare_rth  = h_bare_rth;
g_bare_rph  = h_bare_rph;
g_bare_thth = 1.0 + h_bare_thth;
g_bare_thph = h_bare_thph;
g_bare_phph = 1.0 + h_bare_phph;

g_rr        = g_bare_rr;
g_rth       = rr * g_bare_rth;
g_rph       = rr * sin_theta * g_bare_rph;
g_thth      = rr * rr * g_bare_thth;
g_thph      = rr * rr * sin_theta * g_bare_thph;
g_phph      = rr * rr * sin_theta * sin_theta * g_bare_phph;


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

// //Check if any variable is NaN and report which one
//     if (isnan(A2)) printf("Error: A2 is NaN\n");
//     if (isnan(B2)) printf("Error: B2 is NaN\n");
//     if (isnan(C2)) printf("Error: C2 is NaN\n");
//     if (isnan(sin_theta)) printf("Error:sin is NaN\n");
//     if (isnan(cos_theta)) printf("Error: cos is NaN\n");
//     if (isnan(g_xz)) printf("Error: g_xz is NaN\n");
//     if (isnan(detg)) printf("Error: detg is NaN\n");

// if (detg < 0) {
//     std::cout << "The determinant = " << detg << std::endl;
//     std::cout << " rr = " << rr << std::endl;
//     std::cout << " x = " << xbar << std::endl;
//     std::cout << " y = " << ybar << std::endl;
//     std::cout << " z = " << zbar << std::endl;
// 	std::cout << " g_xx = " << g_xx << std::endl;
//     std::cout << " g_xy = " << g_xy << std::endl;
//     std::cout << " g_xz = " << g_xz << std::endl;
//     std::cout << " g_yy = " << g_yy << std::endl;
//     std::cout << " g_yz = " << g_yz << std::endl;
//     std::cout << " g_zz = " << g_zz<< std::endl;

    // NOTE: assume that var only uses 0-5 instead of VAR::U_SYMGT0-5 to save on
    // memory
    var[0] = g_xx;
    var[1] = g_xy;
    var[2] = g_xz;
    var[3] = g_yy;
    var[4] = g_yz;
    var[5] = g_zz;


    // var[VAR::U_SYMGT0] = g_xx;
    // var[VAR::U_SYMGT1] = g_xy;
    // var[VAR::U_SYMGT2] = g_xz;
    // var[VAR::U_SYMGT3] = g_yy;
    // var[VAR::U_SYMGT4] = g_yz;
    // var[VAR::U_SYMGT5] = g_zz;
}

}  // namespace bssn

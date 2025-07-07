// double mass   = 1.0; 
// double a_spin = 0.1 ; 
// double q_elec = 0.0 ;
double mass   = bssn::BH1.getBHMass();
double a_spin = bssn::BH1.getBHSpin();
double q_elec = bssn::BH1.getBHCharge();
double p_magn = 0.0; 
double rbar ; 
double cyl_rho_bar ; 
double xbar = xx ; 
double ybar = yy ; 
double zbar = zz ; 

double cyl_rho_bar_sqr = xbar*xbar + ybar*ybar ; 
double rbar_sqr = cyl_rho_bar_sqr + zbar*zbar ; 

double cos_theta ;
double sin_theta ; 


if ( fabs(xbar) > 1.0e-4 || fabs(ybar) > 1.0e-4 ) // not on axis  
{ 
    cyl_rho_bar = sqrt( cyl_rho_bar_sqr ) ; 
    rbar = sqrt( rbar_sqr ) ; 
    cos_theta = zbar / rbar ; 
    sin_theta = cyl_rho_bar / rbar ; 
}     
else // on axis  
{
    cyl_rho_bar = 1.0e-10 ; 

    if ( fabs(zbar) > 1.0e-4 ) // on axis but not at the origin  
    { 
        rbar = sqrt( rbar_sqr ) ; 
        cos_theta = zbar / rbar ; 
        sin_theta = 0.0 ; //cyl_rho_bar / rbar ; 
    } 
    else // at the origin  
    {
        rbar = 1.0e-10 ; 
        cos_theta = 0.0 ; 
        sin_theta = 0.0 ; 
    }
}


double r_boylin; 
double f_0 ; 
double rho_bar_sqr ; 
double Delta_bar ; 
double Sigma_bar ; 
double CC ; 
double beta_phi_up ; 


f_0 =  (rbar + 0.5*mass)*(rbar + 0.5*mass) - 0.25*(pow(a_spin,2) + pow(q_elec,2) + pow(p_magn,2) ) ;

rho_bar_sqr = f_0*f_0 + a_spin*a_spin*zbar*zbar ;

Delta_bar = pow(f_0,2) - 2*mass*rbar*f_0 + (a_spin*a_spin+pow(q_elec,2)+pow(p_magn,2))*rbar*rbar ;

Sigma_bar =   pow(f_0*f_0 + a_spin*a_spin*rbar*rbar , 2 ) 
            - a_spin*a_spin*(xbar*xbar + ybar*ybar)*Delta_bar ;

CC = rho_bar_sqr*rho_bar_sqr / Sigma_bar ;

beta_phi_up = - a_spin * pow(rbar,2) * (pow(f_0,2) + pow(a_spin*rbar,2) - Delta_bar ) / Sigma_bar ;

r_boylin = f_0 / rbar ;

var[VAR::U_ALPHA] = sqrt( Delta_bar * rho_bar_sqr / Sigma_bar ) ;
var[VAR::U_CHI] = pow(rbar, 4)/pow(Sigma_bar*rho_bar_sqr, 1.0/3.0) ;
var[VAR::U_K] = 0.0 ;

var[VAR::U_BETA0] = - ybar * beta_phi_up ;    // shift_x  
var[VAR::U_BETA1] =   xbar * beta_phi_up ;    // shift_y   
var[VAR::U_BETA2] = 0.0 ;    // shift_z 

var[VAR::U_GT0] = pow(a_spin,2) * xbar * rho_bar_sqr * pow(CC,-4.0/3.0)/3.0 / pow(Sigma_bar,2) * ( 2 * pow(sin_theta,2) * sqrt( Delta_bar ) * ( rho_bar_sqr * ( f_0 - mass*rbar ) + 2*f_0*rbar*(2*mass*f_0 - rbar*(pow(q_elec,2) + pow(p_magn,2))) ) + 2*pow(cos_theta,2) * (rho_bar_sqr*Delta_bar - 2*Sigma_bar) - 3*rho_bar_sqr*(rho_bar_sqr + 2*mass*f_0*rbar - (pow(q_elec,2)+pow(p_magn,2))*pow(rbar,2)) ) ;

//var[VAR::U_CAP_GT0] = //pow(a_spin,2) * xbar * rho_bar_sqr * pow(CC,-4.0/3.0)/3.0 / pow(Sigma_bar,2) * 
//( 
//2 * pow(sin_theta,2) * sqrt( Delta_bar ) * ( rho_bar_sqr * ( f_0 - mass*rbar ) + 2*f_0*rbar*(2*mass*f_0 - rbar*(pow(q_elec,2) + pow(p_magn,2))) ) 
//+ 2*pow(cos_theta,2) * (rho_bar_sqr*Delta_bar - 2*Sigma_bar) - 3*rho_bar_sqr*(rho_bar_sqr + 2*mass*f_0*rbar - (pow(q_elec,2)+pow(p_magn,2))*pow(rbar,2)) 
//) ;
//;

var[VAR::U_GT1] = pow(a_spin,2) * ybar * rho_bar_sqr * pow(CC,-4.0/3.0)/3.0 / pow(Sigma_bar,2) * ( 2 * pow(sin_theta,2) * sqrt( Delta_bar ) * ( rho_bar_sqr * ( f_0 - mass*rbar ) + 2*f_0*rbar*(2*mass*f_0 - rbar*(pow(q_elec,2) + pow(p_magn,2))) ) + 2*pow(cos_theta,2) * (rho_bar_sqr*Delta_bar - 2 *Sigma_bar ) - 3*rho_bar_sqr*(rho_bar_sqr + 2*mass*f_0*rbar - (pow(q_elec,2)+pow(p_magn,2))*pow(rbar,2)) ) ;

//var[VAR::U_CAP_GT1] = //pow(a_spin,2) * ybar * rho_bar_sqr * pow(CC,-4.0/3.0)/3.0 / pow(Sigma_bar,2) * 
//( 
//2 * pow(sin_theta,2) * sqrt( Delta_bar ) * ( rho_bar_sqr * ( f_0 - mass*rbar ) + 2*f_0*rbar*(2*mass*f_0 - rbar*(pow(q_elec,2) + pow(p_magn,2))) ) +  
//2*pow(cos_theta,2) * (rho_bar_sqr*Delta_bar - 2 *Sigma_bar ) - 3*rho_bar_sqr*(rho_bar_sqr + 2*mass*f_0*rbar - (pow(q_elec,2)+pow(p_magn,2))*pow(rbar,2))   
//) ;
//;

var[VAR::U_GT2] = 2*pow(a_spin,2) * zbar * pow(sin_theta,2) * rho_bar_sqr * pow(CC,-4.0/3.0)/3.0 / pow(Sigma_bar,2) * ( sqrt( Delta_bar ) * ( rho_bar_sqr * ( f_0 - mass*rbar ) + 2*f_0*rbar*(2*mass*f_0 - (pow(q_elec,2)+pow(p_magn,2))*rbar)) - (rho_bar_sqr*Delta_bar - 2 *Sigma_bar ) ) ;

var[VAR::U_B0] = 0.0 ;  // gaugeB0
var[VAR::U_B1] = 0.0 ;  // gaugeB1
var[VAR::U_B2] = 0.0 ;  // gaugeB2
//var[VAR::U_GAUGEB0] = Sigma_bar ;  // gaugeB0
//var[VAR::U_GAUGEB1] = Delta_bar ;  // gaugeB1
//var[VAR::U_GAUGEB2] = CC ;  // gaugeB2

//var[VAR::U_GT00] = pow(CC, -2.0/3.0) * (xbar*xbar*CC + ybar*ybar) / cyl_rho_bar ;
var[VAR::U_SYMGT0] = pow(CC, -2.0/3.0) * (CC - (CC-1.0)*ybar*ybar/cyl_rho_bar) ;
var[VAR::U_SYMGT1] = xbar*ybar * pow(CC, -2.0/3.0) * (CC - 1.0) / cyl_rho_bar ;
var[VAR::U_SYMGT2] = 0.0 ;
//var[VAR::U_GT11] = pow(CC, -2.0/3.0) * (ybar*ybar*CC + xbar*xbar) / cyl_rho_bar ;
var[VAR::U_SYMGT3] = pow(CC, -2.0/3.0) * (CC - (CC-1.0)*xbar*xbar/cyl_rho_bar) ;
var[VAR::U_SYMGT4] = 0.0 ;
var[VAR::U_SYMGT5] = pow(CC, 1.0/3.0) ;

var[VAR::U_SYMAT0] = 2*a_spin*xbar*ybar*rbar/rho_bar_sqr/pow(Sigma_bar*rho_bar_sqr,5.0/6.0) * ( mass*rho_bar_sqr*(pow(f_0,2) + pow(a_spin*rbar,2)) - f_0*(2*mass*f_0 - rbar*(pow(q_elec,2) + pow(p_magn,2)))*(rho_bar_sqr+pow(f_0,2) + pow(a_spin*rbar,2) ) + pow(a_spin*zbar,2)*sqrt(Delta_bar)*(2*mass*f_0 - rbar*(pow(q_elec,2) + pow(p_magn,2))) ) ;

var[VAR::U_SYMAT1] = - a_spin*(xbar*xbar - ybar*ybar)*rbar/rho_bar_sqr/pow(Sigma_bar*rho_bar_sqr,5.0/6.0) * ( mass*rho_bar_sqr*(pow(f_0,2) + pow(a_spin*rbar,2)) - f_0*(2*mass*f_0 - rbar*(pow(q_elec,2) + pow(p_magn,2)))*(rho_bar_sqr+pow(f_0,2) + pow(a_spin*rbar,2) ) + pow(a_spin*zbar,2)*sqrt(Delta_bar)*(2*mass*f_0 - rbar*(pow(q_elec,2) + pow(p_magn,2))) ) ;

var[VAR::U_SYMAT2] = a_spin*ybar*zbar*rbar/rho_bar_sqr/pow(Sigma_bar*rho_bar_sqr,5.0/6.0) * ( mass*rho_bar_sqr*(pow(f_0,2) + pow(a_spin*rbar,2)) - f_0*(2*mass*f_0 - rbar*(pow(q_elec,2)+pow(p_magn,2)))*(rho_bar_sqr+pow(f_0,2) + pow(a_spin*rbar,2) ) - pow(a_spin,2)*(xbar*xbar+ybar*ybar)*sqrt(Delta_bar)*(2*mass*f_0 - rbar*(pow(q_elec,2) + pow(p_magn,2))) ) ;

var[VAR::U_SYMAT3] = - var[VAR::U_SYMAT0];

var[VAR::U_SYMAT4] = - a_spin*xbar*zbar*rbar/rho_bar_sqr/pow(Sigma_bar*rho_bar_sqr,5.0/6.0) * ( mass*rho_bar_sqr*(pow(f_0,2) + pow(a_spin*rbar,2)) - f_0*(2*mass*f_0 - rbar*(pow(q_elec,2)+pow(p_magn,2)))*(rho_bar_sqr+pow(f_0,2) + pow(a_spin*rbar,2) ) - pow(a_spin,2)*(xbar*xbar+ybar*ybar)*sqrt(Delta_bar)*(2*mass*f_0 - rbar*(pow(q_elec,2) + pow(p_magn,2))) ) ;

var[VAR::U_SYMAT5] = 0.0 ;


var[VAR::U_DAMP_PSI] = 0.0 ;  // dampingPsi
var[VAR::U_DAMP_PHI] = 0.0 ;  // dampingPhi
//var[VAR::U_DAMPINGPSI] = var[VAR::U_AT00] ;  // dampingPsi
//var[VAR::U_DAMPINGPHI] = var[VAR::U_AT01] ;  // dampingPhi
//var[VAR::U_DAMPINGPSI] = pow(Sigma_bar,1.0/6.0) ;  // dampingPsi
//var[VAR::U_DAMPINGPHI] = CC-1.0 ;  // dampingPhi

var[VAR::U_PERP_E0] = xbar * pow(rbar,3)/pow(rho_bar_sqr,2)/sqrt(Sigma_bar*rho_bar_sqr) * ( (pow(f_0,2) + pow(a_spin*rbar,2))* ( q_elec*(pow(f_0,2) - pow(a_spin*zbar,2) ) - 2*p_magn*a_spin*f_0*zbar ) - sqrt(Delta_bar)*a_spin*zbar*(2*q_elec*a_spin*f_0*zbar + p_magn*(pow(f_0,2) - pow(a_spin*zbar,2))) ) ;   

var[VAR::U_PERP_E1] = ybar * pow(rbar,3)/pow(rho_bar_sqr,2)/sqrt(Sigma_bar*rho_bar_sqr) * ( (pow(f_0,2) + pow(a_spin*rbar,2))* ( q_elec*(pow(f_0,2) - pow(a_spin*zbar,2) ) - 2*p_magn*a_spin*f_0*zbar ) - sqrt(Delta_bar)*a_spin*zbar*(2*q_elec*a_spin*f_0*zbar + p_magn*(pow(f_0,2) - pow(a_spin*zbar,2))) ) ;   

var[VAR::U_PERP_E2] = pow(rbar,3)/pow(rho_bar_sqr,2)/sqrt(Sigma_bar*rho_bar_sqr) * ( zbar*(pow(f_0,2) + pow(a_spin*rbar,2))* ( q_elec*(pow(f_0,2) - pow(a_spin*zbar,2) ) - 2*p_magn*a_spin*f_0*zbar ) + sqrt(Delta_bar)*a_spin*(xbar*xbar+ybar*ybar)*(2*q_elec*a_spin*f_0*zbar + p_magn*(pow(f_0,2) - pow(a_spin*zbar,2))) ) ;   

var[VAR::U_PERP_B0] = xbar * pow(rbar,3) / pow(rho_bar_sqr,2) / sqrt(Sigma_bar*rho_bar_sqr) * ( ( pow(f_0,2) + pow(a_spin*rbar,2)) * ( 2*q_elec*a_spin*f_0*zbar + p_magn*( pow(f_0,2) - pow(a_spin*zbar,2) ) ) + sqrt(Delta_bar) * a_spin * zbar * ( q_elec*( pow(f_0,2) - pow(a_spin*zbar,2) ) - 2*p_magn*a_spin*f_0*zbar ) ) ;

var[VAR::U_PERP_B1] = ybar * pow(rbar,3)/pow(rho_bar_sqr,2)/sqrt(Sigma_bar*rho_bar_sqr) * ( (pow(f_0,2) + pow(a_spin*rbar,2))* ( 2*q_elec*a_spin*f_0*zbar + p_magn*(pow(f_0,2) - pow(a_spin*zbar,2) ) ) + sqrt(Delta_bar)*a_spin*zbar*( q_elec*(pow(f_0,2) - pow(a_spin*zbar,2)) - 2*p_magn*a_spin*f_0*zbar) ) ;  

var[VAR::U_PERP_B2] = pow(rbar,3)/pow(rho_bar_sqr,2)/sqrt(Sigma_bar*rho_bar_sqr) * ( zbar * (pow(f_0,2) + pow(a_spin*rbar,2))* ( 2*q_elec*a_spin*f_0*zbar + p_magn*(pow(f_0,2) - pow(a_spin*zbar,2) ) ) - sqrt(Delta_bar)*a_spin*(xbar*xbar+ybar*ybar)*( q_elec*(pow(f_0,2) - pow(a_spin*zbar,2)) - 2*p_magn*a_spin*f_0*zbar) ) ;  

// var[VAR::U_PERP_E0] = 0.0;
// var[VAR::U_PERP_E1] = 0.0;
// var[VAR::U_PERP_E2] = 0.0;
// var[VAR::U_PERP_B0] = 0.0; 
// var[VAR::U_PERP_B1] = 0.0;
// var[VAR::U_PERP_B2] = 0.0;


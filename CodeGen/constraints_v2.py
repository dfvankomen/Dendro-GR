####################################################################
#
# Date : Dec.12.2017
# Updated : May.18.2018 
# Python script that generates Psi4 for gravitational waves and
# the momentum and Hamiltonian constraint equations. 
# 
# ewh: updated Jan 2024 to include Riemann^2 (the Kretschmann invariant)
#
#This constraints_v2.py file includes the new diagnostics: Riemann squared, the Pontryagin index and the expansion
####################################################################

#!/usr/bin/env/ python3

import dendro
from sympy import *
dendro.e_ijk  = [(0, 0, 0), (0, 0, 1), (0, 0, 2),
         (0, 1, 0), (0, 1, 1), (0, 1, 2),
         (0, 2, 0), (0, 2, 1), (0, 2, 2),
         (1, 0, 0), (1, 0, 1), (1, 0, 2),
         (1, 1, 0), (1, 1, 1), (1, 1, 2),
         (1, 2, 0), (1, 2, 1), (1, 2, 2),
         (2, 0, 0), (2, 0, 1), (2, 0, 2),
         (2, 1, 0), (2, 1, 1), (2, 1, 2),
         (2, 2, 0), (2, 2, 1), (2, 2, 2)]


###################################################################
# initialize
###################################################################

# Declare variables.
# These include the BSSN variables that we need for the Psi4 
# calculation. 
chi = dendro.scalar("chi","[pp]")
K   = dendro.scalar("K","[pp]")
Gt  = dendro.vec3("Gt","[pp]")
gt  = dendro.sym_3x3("gt","[pp]")
At  = dendro.sym_3x3("At","[pp]")

# Specify the operators needed for computing first and second derivatives
d = dendro.set_first_derivative('grad')    # first argument is direction
d2 = dendro.set_second_derivative('grad2')  # first 2 arguments are directions
ad = dendro.set_advective_derivative('agrad')  # first argument is direction

# Metric related quantities, i.e. the metric and its inverse  
dendro.set_metric(gt)
igt = dendro.get_inverse_metric()

# Christoffels, Ricci, et al  
C1 = dendro.get_first_christoffel()
# Recall that C2 (as C1) is made from the conformally rescaled metric, i.e. it is 
# the "tilded" version of the second Christoffel (or first for C1) symbol.   
C2 = dendro.get_second_christoffel()
#what's this...tried to comment it out and python compilation fails 
C2_spatial = dendro.get_complete_christoffel(chi)
R, Rt, Rphi, CalGt = dendro.compute_ricci(Gt, chi)


###################################################################
# Calculate the tetrad used in the Psi4 calculation
####################################################################

# Define coordinates
x, y, z = symbols('x, y, z')

# Some other values
invsqrt2 = 1/sqrt(2) #0.7071067811865475244
inv_chi = 1/chi

# Define the original spatial vectors in our tetrad 
r_vec = Matrix([[x,y,z]])
theta = Matrix([[x*z,y*z,-(x*x+y*y)]])
phi = Matrix([[-y,x,0.0]])

# We use Gram-Schmidt to make the basis orthonormal. 
# Note that we use the original (not conformally rescaled) metric to define
# the tetrad and correspondingly Psi4.  
gd = gt*inv_chi

# For r_vec
inner_product = 0.0
inner_product = sum([sum([gd[i,j] * r_vec[i] * r_vec[j] for i in dendro.e_i]) for j in dendro.e_i])

r_vec /= sqrt(inner_product)

# For theta
inner_product_1 = 0.0
inner_product_2 = 0.0

inner_product_1 = sum([sum([gd[i,j] * theta[i] * theta[j] for i in dendro.e_i]) for j in dendro.e_i])
inner_product_2 = sum([sum([gd[i,j] * r_vec[i] * theta[j] for i in dendro.e_i]) for j in dendro.e_i])

theta -= inner_product_2 * r_vec
theta /= sqrt(inner_product_1 - inner_product_2 * inner_product_2)

# For phi
inner_product_1 = 0.0
inner_product_2 = 0.0
inner_product_3 = 0.0

inner_product_1 = sum([sum([gd[i,j] * phi[i] * phi[j] for i in dendro.e_i]) for j in dendro.e_i])
inner_product_2 = sum([sum([gd[i,j] * r_vec[i] * phi[j] for i in dendro.e_i]) for j in dendro.e_i])
inner_product_3 = sum([sum([gd[i,j] * theta[i] * phi[j] for i in dendro.e_i]) for j in dendro.e_i])

phi -= inner_product_2 * r_vec + inner_product_3 * theta
phi /= sqrt(inner_product_1 - inner_product_2 * inner_product_2 - inner_product_3 * inner_product_3)

# This completes the tetrad construction. 

###################################################################
# Calculate the Weyl scalar, Psi4, for gravitational wave extraction
###################################################################

# Rename the tetrad quantities for calculating Psi4 
r_np = Matrix([[r_vec[0],r_vec[1],r_vec[2]]])
m_np_real = Matrix([[theta[0],theta[1],theta[2]]])*invsqrt2
m_np_imag = Matrix([[phi[0],phi[1],phi[2]]])*invsqrt2

# Some auxilary variables
# MM and NN are symmetric 2nd rank objects and
# MR and NR are anti-symmetric 2nd rank objects 

MM = Matrix([m_np_real[i]*m_np_real[j] - m_np_imag[i]*m_np_imag[j] for i,j in dendro.e_ij]) 
MM = MM.reshape(3,3)
NN = Matrix([m_np_real[i]*m_np_imag[j] + m_np_real[j]*m_np_imag[i] for i,j in dendro.e_ij])
NN = NN.reshape(3,3)  
MR = Matrix([m_np_real[i]*r_np[j] - m_np_real[j]*r_np[i] for i,j in dendro.e_ij])
MR = MR.reshape(3,3) 
NR = Matrix([m_np_imag[i]*r_np[j] - m_np_imag[j]*r_np[i] for i,j in dendro.e_ij])
NR = NR.reshape(3,3)  

# Additional intermediate variables
#A_vec = Matrix([[sum([At[j,0]*r_np[j] for j in dendro.e_i]), sum([At[j,1]*r_np[j] for j in dendro.e_i]),sum([At[j,2]*r_np[j] for j in dendro.e_i])]])
A_vec = [ sum([At[i,j]*r_np[j] for j in dendro.e_i]) for i in dendro.e_i ] 

Uu = Matrix([sum([m_np_real[k] * (d(j, At[k,i]) + sum([C2[m,k,i] * At[m,j] for m in dendro.e_i])) for k in dendro.e_i]) for i,j in dendro.e_ij])
Uu = Uu.reshape(3,3)
Vv = Matrix([sum([m_np_imag[k] * (d(j, At[k,i]) + sum([C2[m,k,i] * At[m,j] for m in dendro.e_i])) for k in dendro.e_i]) for i,j in dendro.e_ij])
Vv = Vv.reshape(3,3)

r_d_chi = sum([r_np[i] * d(i, chi) for i in dendro.e_i]) 

A_temp = inv_chi * inv_chi * ( sum([A_vec[i] * r_np[i] for i in dendro.e_i]) + K * chi * Rational(1,3) + Rational(1,2) * r_d_chi ) 

m_real_d_chi = sum([m_np_real[i] * d(i, chi) for i in dendro.e_i])  
m_imag_d_chi = sum([m_np_imag[i] * d(i, chi) for i in dendro.e_i]) 

m_real_A_vec = sum([m_np_real[i] * A_vec[i] for i in dendro.e_i]) 
m_imag_A_vec = sum([m_np_imag[i] * A_vec[i] for i in dendro.e_i])  


# Calculate Psi4

psi4_1_real = sum([R[i,i] * MM[i,i] for i in dendro.e_i]) + 2*sum([sum([R[i,j] * MM[i,j] for j in range(i+1,3)]) for i in range(0,2)]) 
psi4_1_imag = sum([R[i,i] * NN[i,i] for i in dendro.e_i]) + 2*sum([sum([R[i,j] * NN[i,j] for j in range(i+1,3)]) for i in range(0,2)]) 

psi4_2_real = A_temp * (sum([At[i,i] * MM[i,i] for i in dendro.e_i]) + 2*sum([sum([At[i,j] * MM[i,j] for j in range(i+1,3)]) for i in range(0,2)]))
psi4_2_imag = A_temp * (sum([At[i,i] * NN[i,i] for i in dendro.e_i]) + 2*sum([sum([At[i,j] * NN[i,j] for j in range(i+1,3)]) for i in range(0,2)]))  

psi4_3_real = inv_chi * sum([sum([MR[i,j]* Uu[i,j] - NR[i,j]*Vv[i,j] for i in dendro.e_i]) for j in dendro.e_i])  
psi4_3_imag = inv_chi * sum([sum([NR[i,j]* Uu[i,j] + MR[i,j]*Vv[i,j] for i in dendro.e_i]) for j in dendro.e_i]) 

#12/31/2020: 0.5 is replaced with rational. 
psi4_4_real = inv_chi * inv_chi * (m_real_A_vec * (m_real_A_vec + Rational(1,2) * m_real_d_chi) - m_imag_A_vec * (m_imag_A_vec + Rational(1,2) * m_imag_d_chi))  
# EWH 13May2024  This is wrong.  Fixed below. 
#psi4_4_imag = inv_chi * inv_chi * (m_real_A_vec * (m_imag_A_vec - Rational(1,2) * m_imag_d_chi ) + m_imag_A_vec * (m_real_A_vec - Rational(1,2) * m_real_d_chi))  
psi4_4_imag = inv_chi * inv_chi * (m_real_A_vec * (m_imag_A_vec + Rational(1,2) * m_imag_d_chi ) + m_imag_A_vec * (m_real_A_vec + Rational(1,2) * m_real_d_chi))  

# Adding previous auxilary Psi4 calculations
# 12/31/2020 : There is a - sign convention issue to match the sign with the LazEv Code. 
#psi4_real =     psi4_1_real + psi4_2_real - psi4_3_real - psi4_4_real
#psi4_img = - ( psi4_1_imag + psi4_2_imag - psi4_3_imag - psi4_4_imag  )
psi4_real = - ( psi4_1_real + psi4_2_real - psi4_3_real - psi4_4_real )
psi4_img =   ( psi4_1_imag + psi4_2_imag - psi4_3_imag - psi4_4_imag )

###################################################################
# Calculate two other Weyl scalars, Psi0 (for ingoing GW radiation) 
# and Psi2 (the Coulomb piece)
###################################################################
# # The calculation of Psi_0 is very similar to that for Psi4. Indeed, 
# # we can write 
# #   $\Psi_0 = \Psi4^{*}(r^i \rightarrow - r^i)$ 
# # and that's it.  So we need to change the sign of all quantities 
# # linear in the vector $r^i$ and swap the sign of the imaginary 
# # part.  (This is the change from outgoing to ingoing radiation.) 
# MR_0 = - MR 
# MR_0 = MR_0.reshape(3,3) 
# NR_0 = - NR 
# NR_0 = NR_0.reshape(3,3)  
# A_vec_0 = -A_vec 
# r_d_chi_0 = - r_d_chi 
# A_temp_0 = inv_chi * inv_chi * ( sum([A_vec[i] * r_np[i] for i in dendro.e_i]) - K * chi * Rational(1,3) + Rational(1,2) * r_d_chi ) 
# m_real_A_vec_0 = - m_real_A_vec 
# m_imag_A_vec_0 = - m_imag_A_vec 
# psi0_1_real = psi4_1_real 
# psi0_1_imag = psi4_1_imag 
# psi0_2_real = A_temp_0 * (sum([At[i,i] * MM[i,i] for i in dendro.e_i]) + 2*sum([sum([At[i,j] * MM[i,j] for j in range(i+1,3)]) for i in range(0,2)]))
# psi0_2_imag = A_temp_0 * (sum([At[i,i] * NN[i,i] for i in dendro.e_i]) + 2*sum([sum([At[i,j] * NN[i,j] for j in range(i+1,3)]) for i in range(0,2)]))  
# psi0_3_real = - psi4_3_real 
# psi0_3_imag = - psi4_3_imag 
# psi0_4_real = inv_chi * inv_chi * (m_real_A_vec_0 * (m_real_A_vec_0 + Rational(1,2) * m_real_d_chi) - m_imag_A_vec_0 * (m_imag_A_vec + Rational(1,2) * m_imag_d_chi))  
# psi0_4_imag = inv_chi * inv_chi * (m_real_A_vec_0 * (m_imag_A_vec_0 + Rational(1,2) * m_imag_d_chi) + m_imag_A_vec_0 * (m_real_A_vec + Rational(1,2) * m_real_d_chi))  

# psi0_real = - ( psi4_1_real + psi4_2_real - psi4_3_real - psi4_4_real )
# psi0_imag = - ( psi4_1_imag + psi4_2_imag - psi4_3_imag - psi4_4_imag )

#

R_rr = sum([ R[i,i] * r_np[i] * r_np[i] for i in dendro.e_i ]) + 2*sum([ sum([ R[i,j] * r_np[i] * r_np[j] for j in range(i+1,3) ]) for i in range(0,2) ]) 
A_temp_2 = inv_chi * sum([A_vec[i] * r_np[i] for i in dendro.e_i]) - Rational(2,3) * K 
At_MM = sum([ At[i,i] * MM[i,i] for i in dendro.e_i ]) + 2*sum([ sum([ At[i,j] * MM[i,j] for j in range (i+1,3) ]) for i in range(0,2) ])  
At_NN = sum([ At[i,i] * NN[i,i] for i in dendro.e_i ]) + 2*sum([ sum([ At[i,j] * NN[i,j] for j in range (i+1,3) ]) for i in range(0,2) ])  



#psi2_real 
#psi2_imag 


###################################################################
# Constraint Equations
###################################################################

# The Hamiltonian constraint
ham = sum(chi*igt[j,k]*R[j,k] for j,k in dendro.e_ij) - dendro.sqr(At) + Rational(2,3)*K**2

# The momentum  constraints 
mom = Matrix([sum([igt[j,k]*(  d(k,At[i,j]) - \
              sum(dendro.C2[m,k,i]*At[j,m] for m in dendro.e_i)) \
                  for j,k in dendro.e_ij]) for i in dendro.e_i]) - \
      Matrix([sum([Gt[j]*At[i,j] for j in dendro.e_i]) for i in dendro.e_i]) -\
      Rational(3,2)*Matrix([ \
            sum([igt[j,k]*At[k,i]*d(j,chi)/chi for j,k in dendro.e_ij])  \
            for i in dendro.e_i]) -\
      Rational(2,3)*Matrix([d(i,K) for i in dendro.e_i])
mom = [item for sublist in mom.tolist() for item in sublist]

# Output for this should be included psi4_real and psi4_img as double precision


########################################################################
# Calculation of the Kretschmann invariant (=Riemann^2) 
#    = [Riemann]^{abcd}.[Riemann]_{abcd} 
########################################################################
# Again, we will use the full (non-conformal) metric here in a few spots. 
# Recall that above, we defined gd = inv_chi*gt. 
# We will also need igd (the inverse metric): 

igd = chi*igt 

# Some of our quadratic terms will involve the 3D spatial Ricci tensor, R_ab, 
# (built out of the non-conformal(!) metric).  In dendro.py, we define the  
# following tensors:  R_ab (full Ricci), Rt_ab (tilde Ricci -- made from the 
# conformal metric, gt) and Rphi_ab (a tensor made from conformal derivatives 
# of the conformal factor).  The key relation is simply that the latter two 
# sum to give the full Ricci:  R_ab = Rt_ab + Rphi_ab.  It is this R_ab 
# that shows up in the calculation of Kretschmann 
#
# We will also need the mixed form for this, or R^a_b.  In particular, we 
# need:  
mixed_Ricci = Matrix([ sum([ R[i,j]*igd[i,k] for i in dendro.e_i ]) for j,k in dendro.e_ij ]) 
mixed_Ricci = mixed_Ricci.reshape(3,3) 

# We need the Ricci scalar (the trace of mixed Ricci) 
tr_Ricci = sum([ mixed_Ricci[i,i] for i in dendro.e_i ])  

# The square of the Ricci tensor:   
Ricci_sqrd = Matrix([ sum([ R[i,j] * mixed_Ricci[j,k] for j in dendro.e_i] ) for i,k in dendro.e_ij ]) 
Ricci_sqrd = Ricci_sqrd.reshape(3,3) 

# The trace of Ricci squared:  
tr_Ricci_sqrd = sum([ Ricci_sqrd[i,i]*igd[i,i] for i in dendro.e_i ]) + 2*sum([ sum([ Ricci_sqrd[i,j]*igd[i,j] for j in range(i+1,3) ]) for i in range(0,2) ])

# We will need several nonlinear terms involving the extrinsic 
# curvature, K_ab.  
# We will first write K_ab in terms of At_ab (tilde A_ab), the conformal 
# factor and K=tr(K_ab)
ext_curv_K = inv_chi * Matrix([ At[i,j] + Rational(1,3) * K * gt[i,j] for i,j in dendro.e_ij ])  
ext_curv_K = ext_curv_K.reshape(3,3) 

# We need the mixed extrinsic curvature, or K^a_b: 
mixed_ext_curv_K = Matrix([ sum([ ext_curv_K[i,j]*igd[i,k] for i in dendro.e_i ]) for j,k in dendro.e_ij ]) 
mixed_ext_curv_K = mixed_ext_curv_K.reshape(3,3) 

# We need the square of the extrinsic curvature:  
ext_curv_K_sqrd = Matrix([ sum([ ext_curv_K[i,j] * mixed_ext_curv_K[j,k] for j in dendro.e_i ]) for i,k in dendro.e_ij ]) 
ext_curv_K_sqrd = ext_curv_K_sqrd.reshape(3,3) 

# We actually need the mixed version of K squared (K^k_i K^i_j) 
mixed_ext_curv_K_sqrd = Matrix([ sum([ ext_curv_K_sqrd[i,j]*igd[i,k] for i in dendro.e_i ]) for j,k in dendro.e_ij ]) 
mixed_ext_curv_K_sqrd = mixed_ext_curv_K_sqrd.reshape(3,3) 

# And the cube of the extrinsic curvature:  
ext_curv_K_cubed = Matrix([ sum([ ext_curv_K_sqrd[i,j] * mixed_ext_curv_K[j,k] for j in dendro.e_i ]) for i,k in dendro.e_ij ]) 
ext_curv_K_cubed = ext_curv_K_cubed.reshape(3,3) 

# And the extrinsic curvature to the fourth: 
ext_curv_K_foured = Matrix([ sum([ ext_curv_K_cubed[i,j] * mixed_ext_curv_K[j,k] for j in dendro.e_i ]) for i,k in dendro.e_ij ]) 
ext_curv_K_foured = ext_curv_K_foured.reshape(3,3) 

# As well as the trace of the square, cube and fourth powers:
tr_ext_curv_K_sqrd   = sum([ ext_curv_K_sqrd[i,i]  *igd[i,i] for i in dendro.e_i ]) + 2*sum([ sum([ ext_curv_K_sqrd[i,j]  *igd[i,j] for j in range(i+1,3) ]) for i in range(0,2) ])
tr_ext_curv_K_cubed  = sum([ ext_curv_K_cubed[i,i] *igd[i,i] for i in dendro.e_i ]) + 2*sum([ sum([ ext_curv_K_cubed[i,j] *igd[i,j] for j in range(i+1,3) ]) for i in range(0,2) ])
tr_ext_curv_K_foured = sum([ ext_curv_K_foured[i,i]*igd[i,i] for i in dendro.e_i ]) + 2*sum([ sum([ ext_curv_K_foured[i,j]*igd[i,j] for j in range(i+1,3) ]) for i in range(0,2) ])

# We define a few cross terms between the extrinsic and Ricci curvatures 
# and make them temporary variables  
tmp_1 = Matrix([ sum([ mixed_Ricci[i,j] * ( K * mixed_ext_curv_K[j,k] - mixed_ext_curv_K_sqrd[j,k] ) for j in dendro.e_i ]) for i,k in dendro.e_ij ]) 
tmp_1 = tmp_1.reshape(3,3) 

tr_tmp_1 = sum([ tmp_1[i,i] for i in dendro.e_i ]) 

# Matter pieces with the 3-projected stress-tensor, S_ab, its trace, S, and 
# the energy density, rho.  (We will assume units such that G=1.)  
#
# For the moment, we seem to be using this calculation for gravitational wave 
# (vacuum, Teukolsky) evolutions.  If and when we include matter, these terms
# will need to be included in the full calculation.  To that end, be sure 
# to include YOURFAVORITE* pieces below for the matter terms.  
#
#eight_pi_G = 8*pi  
#proj_T = Matrix([ YOURFAVORITEMATTERSOURCE[i,j] for i,j in dendro.e_ij ])  
#proj_T = proj_T.reshape(3,3) 
#rho = YOURFAVORITEENERGYDENSITY 
#mixed_proj_T = Matrix([ sum([ proj_T[i,j]*igd[i,k] for i in dendro.e_i ]) for j,k in dendro.e_ij ]) 
#mixed_proj_T = mixed_proj_T.reshape(3,3) 
#tr_proj_T = sum([ mixed_proj_T[i,i] for i in dendro.e_i ]) 
#proj_T_sqrd = Matrix([ sum([ proj_T[i,j] * mixed_proj_T[j,k] for j in dendro.e_i] ) for i,k in dendro.e_ij ]) 
#proj_T_sqrd = proj_T_sqrd.reshape(3,3) 
#tr_proj_T_sqrd = sum([ proj_T_sqrd[i,i]*igd[i,i] for i in dendro.e_i ]) + 2*sum([ sum([ proj_T_sqrd[i,j]*igd[i,j] for j in range(i+1,3) ]) for i in range(0,2) ])
# note that the K K^j_k - K_ji K^i_k used here is used other places; optimize? 
#tmp_2 = Matrix([ sum([ mixed_proj_T[i,j] * ( mixed_Ricci[j,k] + K * mixed_ext_curv_K[j,k] - mixed_ext_curv_K_sqrd[j,k] ) for j in dendro.e_i ]) for i,k in dendro.e_ij ]) 
#tmp_2 = tmp_2.reshape(3,3) 
#tr_tmp_2 = sum([ tmp_2[i,i] for i in dendro.e_i ]) 
# We give the first part of Riemann squared.  (This one is for matter spacetimes.) 
#Riemann_sqrd_part1=   8*tr_Ricci_sqrd 
#                     - tr_Ricci*tr_Ricci 
#                     + 16*tr_tmp_1 
#                     - 2*tr_Ricci*(K*K - tr_ext_curv_K_sqrd) 
#                     + 2*tr_ext_curv_K_sqrd*tr_ext_curv_K_sqrd 
#                     + 2*tr_ext_curv_K_foured 
#                     + 4*K*(K*tr_ext_curv_K_sqrd - 2*tr_ext_curv_K_cubed) 
#                     + 4*eight_pi_G*( - 2*tr_tmp_2 + (tr_proj_T-rho)*( tr_Ricci + K*K - tr_ext_curv_K_sqrd) + eight_pi_G *( tr_proj_T_sqrd - Rational(1,4)*(tr_proj_T-rho)*(tr_proj_T+3*rho) ) )  

# We give the first part of Riemann squared.  (This one is for vacuum spacetimes.) 
riemann_sqrd_part1 =   8*tr_Ricci_sqrd - tr_Ricci*tr_Ricci + 16*tr_tmp_1 - 2*tr_Ricci*(K*K - tr_ext_curv_K_sqrd) + 2*tr_ext_curv_K_sqrd*tr_ext_curv_K_sqrd + 2*tr_ext_curv_K_foured + 4*K*(K*tr_ext_curv_K_sqrd - 2*tr_ext_curv_K_cubed)    


# Now to include terms with covariant derivatives of the extrinsic curvature 

#mixed_At = dendro.up_down(At) 

mixed_At = Matrix([ sum([ At[i,j]*igt[i,k] for i in dendro.e_i ]) for j,k in dendro.e_ij ]) 
mixed_At = mixed_At.reshape(3,3) 

At_sqrd = Matrix([ sum([ At[i,j] * mixed_At[j,k] for j in dendro.e_i ]) for i,k in dendro.e_ij ]) 
At_sqrd = At_sqrd.reshape(3,3) 

tr_At_sqrd = sum([ At_sqrd[i,i]*igt[i,i] for i in dendro.e_i ]) + 2*sum([ sum([ At_sqrd[i,j]*igt[i,j] for j in range(i+1,3) ]) for i in range(0,2) ]) 

At_dchi = [ sum([ mixed_At[i,j] * d(i,chi) for i in dendro.e_i ]) for j in dendro.e_i ] 

dchi_sqrd = sum([ d(i,chi)*d(i,chi)*igt[i,i] for i in dendro.e_i ]) + 2*sum([ sum([ d(i,chi)*d(j,chi)*igt[i,j] for j in range(i+1,3) ]) for i in range(0,2) ]) 

dK_sqrd = sum([ d(i,K)*d(i,K)*igt[i,i] for i in dendro.e_i ]) + 2*sum([ sum([ d(i,K)*d(j,K)*igt[i,j] for j in range(i+1,3) ]) for i in range(0,2) ]) 

A_dK_dchi = sum([ sum([ igt[i,j]*At_dchi[i] for i in dendro.e_i ]) * d(j,K) for j in dendro.e_i ])  

# The following is the antisymmetrized (conformally) covariant derivative of "A-tilde" (what we call At_ij), which could be written 
#  ${\tilde D}_{[i} {\tilde A}_{j]k}$ .    
tmp_3 = 0.5 * Array([ d(i,At[j,k]) - d(j,At[i,k]) - sum([ (C2[m,k,i]*At[j,m] - C2[m,k,j]*At[i,m]) for m in dendro.e_i ]) for i,j,k in dendro.e_ijk ]) 
tmp_3 = tmp_3.reshape(3,3,3) 

# The following is basically tmp_3 with the first two indices raised (by the conformal metric, of course):  
tmp_4 = Array([ sum([ sum([ tmp_3[i,j,k] * igt[i,l] for i in dendro.e_i ]) * igt[j,m] for j in dendro.e_i ]) for l,m,k in dendro.e_ijk ]) 
tmp_4 = tmp_4.reshape(3,3,3) 

tmp_5 = Rational(1,2)*inv_chi*Array([ sum([ sum([ (At[k,i]*d(j,chi) - At[k,j]*d(i,chi) )*igt[i,l] for i in dendro.e_i ]) * igt[j,m] for j in dendro.e_i ]) for l,m,k in dendro.e_ijk ])     
tmp_5 = tmp_5.reshape(3,3,3) 

# The following is the divergence of A-tilde:  (TeX-writable as ${\tilde D}_i {\tilde A}_j{}^i$)  
div_At = [ sum([ ( d(i,At[j,k]) - sum([ C2[m,i,j]*At[m,k] + C2[m,i,k]*At[j,m] for m in dendro.e_i ]) ) * igt[i,k] for i,k in dendro.e_ij ]) for j in dendro.e_i ]

tmp_6 = Matrix([ sum([ sum([ tmp_3[i,j,k] * ( tmp_4[i,j,l] + tmp_5[i,j,l] ) for i in dendro.e_i ]) for j in dendro.e_i ]) for k,l in dendro.e_ij ])  
tmp_6 = tmp_6.reshape(3,3) 

tr_tmp_6 = sum([ sum([ tmp_6[i,j]*igt[i,j] for i in dendro.e_i ]) for j in dendro.e_i ]) 

#riemann_sqrd_part2 = - 16*chi*(   tr_tmp_6 
#                                + Rational(1,9)*dK_sqrd 
#                                + Rational(1,8)*tr_At_sqrd*dchi_sqrd/(chi*chi) 
#                                - sum([ sum([ (   div_At[i]*(   Rational(1,2)*inv_chi*At_dchi[j] 
#                                                              + Rational(1,3)*d(j,K) ) 
#                                                - Rational(1,2)*inv_chi*At_dchi[i]*(   d(j,K) 
#                                                                                     + Rational(3,4)*inv_chi*At_dchi[j] 
#                                                                                   ) * igt[i,j] for i in dendro.e_i 
#                                              ) 
#                                            ]) for j in dendro.e_j 
#                                     ]) 
#                              )   

riemann_sqrd_part2 = - 16*chi*(   tr_tmp_6 + Rational(1,9)*dK_sqrd + Rational(1,8)*tr_At_sqrd*dchi_sqrd/(chi*chi) - sum([ sum([ ( div_At[i]*( Rational(1,2)*inv_chi*At_dchi[j] + Rational(1,3)*d(j,K) ) - Rational(1,2)*inv_chi*At_dchi[i]*( d(j,K) + Rational(3,4)*inv_chi*At_dchi[j] ) * igt[i,j] ) for i in dendro.e_i ]) for j in dendro.e_i ]) ) 

riemann_sqrd = riemann_sqrd_part1 + riemann_sqrd_part2  


###################################################################
# Calculation of the Pontryagin index 
#  = [left dual Riemann]^{abcd}.[Riemann]_{abcd}  

levi_civita = [
                 [  
                    [0,  0, 0], 
                    [0,  0, 1], 
                    [0, -1, 0]  
                                ],
                 [  
                    [0, 0, -1], 
                    [0, 0,  0], 
                    [1, 0,  0]  
                                ],
                 [ 
                    [ 0, 1, 0], 
                    [-1, 0, 0], 
                    [ 0, 0, 0]
                                ],
              ]
 
four_pi_G = 4*pi 

#mixed_Ricci = Matrix([ sum([ R[i,j]*igd[i,k] for i in dendro.e_i ]) for j,k in dendro.e_ij ]) 
#mixed_Ricci = mixed_Ricci.reshape(3,3) 

# The first is for matter problems.  
#tmp_pontryagin_1 = Matrix([ mixed_Ricci[i,j] - four_pi_G * mixed_proj_T[i,j] for i,j in dendro.e_ij ])   
tmp_pontryagin_1 = mixed_Ricci
tmp_pontryagin_1 = tmp_pontryagin_1.reshape(3,3) 

#cov_deriv_At[i][j][b] - 2 * sum(At[i,b]*d(phi,j)   

# The next term is just the (conformal) covariant derivative of A-tilde, but with a Christoffel term dropped as it vanishes when contracted against levi-civita as we do at the end.   
tmp_pontryagin_2 = Array([ d(i,At[j,k]) - sum([ C2[m,i,k] * At[j,m] for m in dendro.e_i ]) for i,j,k in dendro.e_ijk ])  
tmp_pontryagin_2 = tmp_pontryagin_2.reshape(3,3,3) 

#tmp_3 = 0.5 * Array([ d(i,At[j,k]) - d(j,At[i,k]) - sum([ (C2[m,k,i]*At[j,m] - C2[m,k,j]*At[i,m]) for m in dendro.e_i ]) for i,j,k in dendro.e_ijk ]) 
tmp_3 = 0.5 * Array([ d(i,At[j,k]) - d(j,At[i,k]) - sum([ (C2[m,k,i]*At[j,m] - C2[m,k,j]*At[i,m]) for m in dendro.e_i ]) for i,j,k in dendro.e_ijk ]) 
tmp_3 = tmp_3.reshape(3,3,3) 

#tmp_pontryagin_3  

tmp_pontryagin_3 = Array([ tmp_pontryagin_2[i,j,k] + Rational(1,2) * At[i,k]*d(chi,j)/chi for i,j,k in dendro.e_ijk ])  
tmp_pontryagin_3 = tmp_pontryagin_3.reshape(3,3,3)  

tmp_pontryagin_31 = Array([ 2*sum([ tmp_pontryagin_2[i,j,b] * tmp_pontryagin_1[b,k] for b in dendro.e_i ]) for i,j,k in dendro.e_ijk ])  
tmp_pontryagin_31 = tmp_pontryagin_31.reshape(3,3,3)  

tmp_pontryagin_4 = Array([ sum([ tmp_pontryagin_2[i,j,b] * ( Rational(2,3) * K * mixed_At[b,k] - sum([ mixed_At[b,e] * mixed_At[e,k] for e in dendro.e_i ]) ) for b in dendro.e_i ]) for i,j,k in dendro.e_ijk ]) 
tmp_pontryagin_4 = tmp_pontryagin_4.reshape(3,3,3) 

tmp_pontryagin_5 = Array([ sum([ sum([ tmp_pontryagin_2[a,b,k] * mixed_At[a,i] for a in dendro.e_i ]) * mixed_At[b,j] for b in dendro.e_i ]) for i,j,k in dendro.e_ijk ]) 
tmp_pontryagin_5 = tmp_pontryagin_5.reshape(3,3,3) 

pontryagin = 8 * sqrt(chi) * sum([ sum([ sum([ levi_civita[i][j][k] * 
( tmp_pontryagin_31[i,j,k] + tmp_pontryagin_4[i,j,k] - tmp_pontryagin_5[i,j,k] ) 
for i in dendro.e_i ]) for j in dendro.e_i ]) for k in dendro.e_i ])   
#pontryagin = 8 * sqrt(chi) * sum([ sum([ sum([ levi_civita[i][j][k] * tmp_pontryagin_4[i,j,k] for i in dendro.e_i ]) for j in dendro.e_i ]) for k in dendro.e_i ])   


# full Christoffel symbols w.r.t g_{ij}
C3 = dendro.get_complete_christoffel(chi)


Kij = (1 / chi) * (At + Rational(1, 3) * gt * K)

expansion = sum([ d(i,r_vec[i]) for i in dendro.e_i])+ sum([sum([C3[i,i,j] for i in dendro.e_i])*r_vec[j]for j in dendro.e_i]) - K + sum([ Kij[i,i]* r_vec[i]*r_vec[i] for i in dendro.e_i]) + 2* sum([sum([Kij[i,j]*r_vec[i]*r_vec[j] for j in range(i+1,3) ]) for i in range(0,2) ])


###################################################################
# generate code
###################################################################
# uncomment to terminal code gen
outs = [psi4_real, psi4_img, ham, mom, riemann_sqrd, pontryagin]
#outs = [psi4_real, psi4_img, ham, mom, riemann_sqrd]
#outs = [psi4_real, psi4_img, ham, mom]
vnames = ['psi4_real', 'psi4_img', 'ham', 'mom', 'riem_sqrd', 'pontryagin']
#vnames = ['psi4_real', 'psi4_img', 'ham', 'mom', 'riem_sqrd']
#vnames = ['psi4_real', 'psi4_img', 'ham', 'mom']
dendro.generate_cpu(outs, vnames, '[pp]')

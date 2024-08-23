####################################################################
#
# Date : Dec.12.2017
# Updated : May.18.2018 
# Python script that generates Psi4 for gravitational waves and
# the momentum and Hamiltonian constraint equations. 
# 
####################################################################

#!/usr/bin/env/ python3

import dendrosym
import sympy as sym
from sympy import *

###################################################################
# initialize
###################################################################

# DEFINE: dendro config class to use for generating code
dendroConfigs = dendrosym.NRConfig("emda")

# the indexing used for the dendro configs
idx_str = "[pp]"

# save the index string
dendroConfigs.set_idx_str(idx_str)
def bssn_constraint_func():
      # Declare variables.
      # These include the BSSN variables that we need for the Psi4 
      # calculation. 
      ham_constraint = dendrosym.dtypes.scalar("ham" + idx_str)
      mom_constraint = dendrosym.dtypes.vec3("mom" + idx_str)
      psi4_real_constraint = dendrosym.dtypes.scalar("psi4_real" + idx_str)
      psi4_imag_constraint = dendrosym.dtypes.scalar("psi4_imag" + idx_str)

      chi = dendrosym.dtypes.scalar("chi","[pp]")
      K   = dendrosym.dtypes.scalar("K","[pp]")
      Gt  = dendrosym.dtypes.vec3("Gt","[pp]")
      gt  = dendrosym.dtypes.sym_3x3("gt","[pp]")
      At  = dendrosym.dtypes.sym_3x3("At","[pp]")

      dendroConfigs.add_constraint_variables(
      [
            ham_constraint,
            mom_constraint,
            psi4_real_constraint,
            psi4_imag_constraint,
      ]
      )

      d_, d2_ = dendrosym.derivs.get_derivs()

      # advective derivate, first argument is direction
      ad_ = dendrosym.nr.set_advective_derivative(d_)
      # and then we set the kreiss oliger dissipation
      kod_ = dendrosym.nr.set_kreiss_oliger_dissipation("kograd")
      # == END DERIVATIVES ==

      # Lie derivative weight
      weight = -sym.Rational(2, 3)
      weight_Gt = sym.Rational(2, 3)

      dendroConfigs.set_metric(gt)
      # and then we get the inverse (conformal) metric
      igt = dendrosym.nr.get_inverse_metric()
      inv_metric = igt
      # as well as the two Christoffel symbols built from the conformal metric
      C1 = dendrosym.nr.get_first_christoffel()
      C2 = dendrosym.nr.get_second_christoffel()
      # and the third Christoffel symbol
      C3 = dendrosym.nr.get_complete_christoffel(chi)

      C2_spatial = C3  # the spatial C2 is the complete christoffel
      R, Rt, Rphi, CalGt = dendrosym.nr.compute_ricci(Gt, chi)

      ###################################################################
      # Calculate the tetrad used in the Psi4 calculation
      ####################################################################

      # Define coordinates
      x, y, z = symbols('x, y, z')

      # Some other values
      invsqrt2 = 1/sqrt(2) #0.7071067811865475244
      inv_chi = 1/chi

      # Define the original spatial vectors in our tetrad 
      r_vec = sym.Matrix([[x,y,z]])
      theta = sym.Matrix([[x*z,y*z,-(x*x+y*y)]])
      phi = sym.Matrix([[-y,x,0.0]])

      # We use Gram-Schmidt to make the basis orthonormal. 
      # Note that we use the original (not conformally rescaled) metric to define
      # the tetrad and correspondingly Psi4.  
      gd = gt*inv_chi

      # For r_vec
      inner_product = 0.0
      inner_product = sum([sum([gd[i,j] * r_vec[i] * r_vec[j] for i in dendrosym.nr.e_i]) for j in dendrosym.nr.e_i])

      r_vec /= sqrt(inner_product)

      # For theta
      inner_product_1 = 0.0
      inner_product_2 = 0.0

      inner_product_1 = sum([sum([gd[i,j] * theta[i] * theta[j] for i in dendrosym.nr.e_i]) for j in dendrosym.nr.e_i])
      inner_product_2 = sum([sum([gd[i,j] * r_vec[i] * theta[j] for i in dendrosym.nr.e_i]) for j in dendrosym.nr.e_i])

      theta -= inner_product_2 * r_vec
      theta /= sqrt(inner_product_1 - inner_product_2 * inner_product_2)

      # For phi
      inner_product_1 = 0.0
      inner_product_2 = 0.0
      inner_product_3 = 0.0

      inner_product_1 = sum([sum([gd[i,j] * phi[i] * phi[j] for i in dendrosym.nr.e_i]) for j in dendrosym.nr.e_i])
      inner_product_2 = sum([sum([gd[i,j] * r_vec[i] * phi[j] for i in dendrosym.nr.e_i]) for j in dendrosym.nr.e_i])
      inner_product_3 = sum([sum([gd[i,j] * theta[i] * phi[j] for i in dendrosym.nr.e_i]) for j in dendrosym.nr.e_i])

      phi -= inner_product_2 * r_vec + inner_product_3 * theta
      phi /= sqrt(inner_product_1 - inner_product_2 * inner_product_2 - inner_product_3 * inner_product_3)

      # This completes the tetrad construction. 

      ###################################################################
      # Calculate the Weyl scalar, Psi4, for graviational wave extraction
      ###################################################################

      # Rename the tetrad quantities for calculating Psi4 
      r_np = sym.Matrix([[r_vec[0],r_vec[1],r_vec[2]]])
      m_np_real = sym.Matrix([[theta[0],theta[1],theta[2]]])*invsqrt2
      m_np_img = sym.Matrix([[phi[0],phi[1],phi[2]]])*invsqrt2

      # Some auxilary variables
      # MM and NN are symmetric 2nd rank objects and
      # MR and NR are anti-symmetric 2nd rank objects 

      MM = sym.Matrix([m_np_real[i]*m_np_real[j] - m_np_img[i]*m_np_img[j] for i,j in dendrosym.nr.e_ij]) 
      MM = MM.reshape(3,3)
      NN = sym.Matrix([m_np_real[i]*m_np_img[j] + m_np_real[j]*m_np_img[i] for i,j in dendrosym.nr.e_ij])
      NN = NN.reshape(3,3)  
      MR = sym.Matrix([m_np_real[i]*r_np[j] - m_np_real[j]*r_np[i] for i,j in dendrosym.nr.e_ij])
      MR = MR.reshape(3,3) 
      NR = sym.Matrix([m_np_img[i]*r_np[j] - m_np_img[j]*r_np[i] for i,j in dendrosym.nr.e_ij])
      NR = NR.reshape(3,3)  

      # Additional intermediate variables
      #A_vec = sym.Matrix([[sum([At[j,0]*r_np[j] for j in dendrosym.nr.e_i]), sum([At[j,1]*r_np[j] for j in dendrosym.nr.e_i]),sum([At[j,2]*r_np[j] for j in dendrosym.nr.e_i])]])
      A_vec = [ sum([At[i,j]*r_np[j] for j in dendrosym.nr.e_i]) for i in dendrosym.nr.e_i ] 

      Uu = sym.Matrix([sum([m_np_real[k] * (d_(j, At[k,i]) + sum([C2[m,k,i] * At[m,j] for m in dendrosym.nr.e_i])) for k in dendrosym.nr.e_i]) for i,j in dendrosym.nr.e_ij])
      Uu = Uu.reshape(3,3)
      Vv = sym.Matrix([sum([m_np_img[k] * (d_(j, At[k,i]) + sum([C2[m,k,i] * At[m,j] for m in dendrosym.nr.e_i])) for k in dendrosym.nr.e_i]) for i,j in dendrosym.nr.e_ij])
      Vv = Vv.reshape(3,3)

      r_d_chi = sum([r_np[i] * d_(i, chi) for i in dendrosym.nr.e_i]) 

      A_temp = inv_chi * inv_chi * ( sum([A_vec[i] * r_np[i] for i in dendrosym.nr.e_i]) + K * chi * Rational(1,3) + Rational(1,2) * r_d_chi ) 

      m_real_d_chi = sum([m_np_real[i] * d_(i, chi) for i in dendrosym.nr.e_i])  
      m_img_d_chi  = sum([m_np_img [i] * d_(i, chi) for i in dendrosym.nr.e_i]) 

      m_real_A_vec = sum([m_np_real[i] * A_vec[i] for i in dendrosym.nr.e_i]) 
      m_img_A_vec  = sum([m_np_img [i] * A_vec[i] for i in dendrosym.nr.e_i])  


      # Calculate Psi4

      psi4_1_real = sum([R[i,i] * MM[i,i] for i in dendrosym.nr.e_i]) + 2*sum([sum([R[i,j] * MM[i,j] for j in range(i+1,3)]) for i in range(0,2)]) 
      psi4_1_img  = sum([R[i,i] * NN[i,i] for i in dendrosym.nr.e_i]) + 2*sum([sum([R[i,j] * NN[i,j] for j in range(i+1,3)]) for i in range(0,2)]) 

      psi4_2_real = A_temp * (sum([At[i,i] * MM[i,i] for i in dendrosym.nr.e_i]) + 2*sum([sum([At[i,j] * MM[i,j] for j in range(i+1,3)]) for i in range(0,2)]))
      psi4_2_img = A_temp * (sum([At[i,i] * NN[i,i] for i in dendrosym.nr.e_i]) + 2*sum([sum([At[i,j] * NN[i,j] for j in range(i+1,3)]) for i in range(0,2)]))  

      psi4_3_real = inv_chi * sum([sum([MR[i,j]* Uu[i,j] - NR[i,j]*Vv[i,j] for i in dendrosym.nr.e_i]) for j in dendrosym.nr.e_i])  
      psi4_3_img = inv_chi * sum([sum([NR[i,j]* Uu[i,j] + MR[i,j]*Vv[i,j] for i in dendrosym.nr.e_i]) for j in dendrosym.nr.e_i]) 

      #12/31/2020: 0.5 is replaced with rational. 
      psi4_4_real = inv_chi * inv_chi * (m_real_A_vec * (m_real_A_vec + Rational(1,2) * m_real_d_chi) - m_img_A_vec * (m_img_A_vec + Rational(1,2) * m_img_d_chi))  
      psi4_4_img = inv_chi * inv_chi * (m_real_A_vec * (m_img_A_vec + Rational(1,2) * m_img_d_chi ) + m_img_A_vec * (m_real_A_vec + Rational(1,2) * m_real_d_chi))  
      # this is potentially wrong but is kept here for comparison purposes. 
      #psi4_4_img = inv_chi * inv_chi * (m_real_A_vec * (m_img_A_vec - Rational(1,2) * m_img_d_chi ) + m_img_A_vec * (m_real_A_vec - Rational(1,2) * m_real_d_chi)) 
      # Adding previous auxilary Psi4 calculations
      # 12/31/2020 : There is a - sign convention issue to match the sign with the LazEv Code. 
      #psi4_real =     psi4_1_real + psi4_2_real - psi4_3_real - psi4_4_real
      #psi4_img  = - ( psi4_1_img  + psi4_2_img  - psi4_3_img  - psi4_4_img  )
      psi4_real =  -(psi4_1_real + psi4_2_real - psi4_3_real - psi4_4_real)
      psi4_img  =  ( psi4_1_img  + psi4_2_img  - psi4_3_img  - psi4_4_img)
      ###################################################################
      # Constraint Equations
      ###################################################################

      # The Hamiltonian constraint
      ham = sum(chi*igt[j,k]*R[j,k] for j,k in dendrosym.nr.e_ij) - dendro.sqr(At) + Rational(2,3)*K**2

      # The momentum  constraints 
      mom = sym.Matrix([sum([igt[j,k]*(  d_(k,At[i,j]) - \
                  sum(dendro.C2[m,k,i]*At[j,m] for m in dendrosym.nr.e_i)) \
                  for j,k in dendrosym.nr.e_ij]) for i in dendrosym.nr.e_i]) - \
      sym.Matrix([sum([Gt[j]*At[i,j] for j in dendrosym.nr.e_i]) for i in dendrosym.nr.e_i]) -\
      Rational(3,2)*sym.Matrix([ \
            sum([igt[j,k]*At[k,i]*d_(j,chi)/chi for j,k in dendrosym.nr.e_ij])  \
            for i in dendrosym.nr.e_i]) -\
      Rational(2,3)*sym.Matrix([d_(i,K) for i in dendrosym.nr.e_i])
      mom = [item for sublist in mom.tolist() for item in sublist]

      # Output for this should be included psi4_real and psi4_img as double precision  
      ###################################################################
      # generate code
      ###################################################################
      #uncomment to terminal code gen
      rhs_list = [psi4_real, psi4_img, ham, mom]
      var_list = [
            psi4_real_constraint,
            psi4_imag_constraint,
            ham_constraint,
            mom_constraint,
      ]

      return rhs_list, var_list


dendroConfigs.set_rhs_equation_function("constraint", bssn_constraint_func)


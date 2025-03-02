#ifndef RHS_H
#define RHS_H

#include <time.h>

#include <cmath>
#include <iostream>

#include "derivs.h"
#include "mathUtils.h"
#include "nlsmUtils.h"
#include "parameters.h"

void nlsmRhs(double **uzipVarsRHS, const double **uZipVars,
             const unsigned int &offset, const double *ptmin,
             const double *ptmax, const unsigned int *sz,
             const unsigned int &bflag);

void nlsm_bcs(double *f_rhs, const double *f, const double *dxf,
              const double *dyf, const double *dzf, const double *pmin,
              const double *pmax, const double f_falloff,
              const double f_asymptotic, const unsigned int *sz,
              const unsigned int &bflag);

void fake_initial_data(double x, double y, double z, double *u);

// void nlsmRhs_blkwise();

#endif

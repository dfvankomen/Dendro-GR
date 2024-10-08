#ifndef RHS_H
#define RHS_H

#include <time.h>

#include <cmath>
#include <iostream>

#include "derivs.h"
#include "mathUtils.h"
#include "parameters.h"
#include "profile_params.h"

#ifdef EM1_ENABLE_CUDA
#include "params_cu.h"
#include "profile_gpu.h"
#include "rhs_cuda.cuh"
#endif

/**@brief computes complete RHS iteratiing over all the blocks.
 * @param[out] unzipVarsRHS: unzipped variables computed RHS
 * @param[in]  unzipVars: unzipped variables.
 * @param[in]  blkList: block list.
 * @param[in]  numBlocks: number of blocks.
 */
void em1rhs(double **uzipVarsRHS, const double **uZipVars,
            const ot::Block *blkList, unsigned int numBlocks);

void em1rhs(double **uzipVarsRHS, const double **uZipVars,
            const unsigned int &offset, const double *ptmin,
            const double *ptmax, const unsigned int *sz,
            const unsigned int &bflag);

void em1_bcs(double *f_rhs, const double *f, const double *dxf,
             const double *dyf, const double *dzf, const double *pmin,
             const double *pmax, const double f_falloff,
             const double f_asymptotic, const unsigned int *sz,
             const unsigned int &bflag);

void em1_bcs_exact_A0(double *f_rhs, const double *f, const double *dxf,
                      const double *dyf, const double *dzf, const double *pmin,
                      const double *pmax, const double f_falloff,
                      const double f_asymptotic, const unsigned int *sz,
                      const unsigned int &bflag);

void em1_bcs_exact_A1(double *f_rhs, const double *f, const double *dxf,
                      const double *dyf, const double *dzf, const double *pmin,
                      const double *pmax, const double f_falloff,
                      const double f_asymptotic, const unsigned int *sz,
                      const unsigned int &bflag);

void em1_bcs_exact_E0(double *f_rhs, const double *f, const double *dxf,
                      const double *dyf, const double *dzf, const double *pmin,
                      const double *pmax, const double f_falloff,
                      const double f_asymptotic, const unsigned int *sz,
                      const unsigned int &bflag);

void em1_bcs_exact_E1(double *f_rhs, const double *f, const double *dxf,
                      const double *dyf, const double *dzf, const double *pmin,
                      const double *pmax, const double f_falloff,
                      const double f_asymptotic, const unsigned int *sz,
                      const unsigned int &bflag);

#endif

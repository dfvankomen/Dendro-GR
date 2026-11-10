//
// Created by milinda on 10/21/17.
/**
 *@author Milinda Fernando
 *School of Computing, University of Utah
 *@brief
 */
//

#include "profile_params.h"
#include "dendroProfileParams.h"

namespace bssn {
namespace timer {

profiler_t total_runtime;

profiler_t t_f2o;
profiler_t t_cons;
profiler_t t_cons_unzip;
profiler_t t_cons_kernel;
profiler_t t_cons_zipex;
profiler_t t_cons_deriv;
profiler_t t_cons_pts;
profiler_t t_bal;
profiler_t t_mesh;

profiler_t t_rkSolve;

profiler_t t_ghostEx_sync;

profiler_t t_unzip_sync;

profiler_t t_unzip_async;

profiler_t t_deriv;
profiler_t t_rhs;
profiler_t t_rhs_ko;
#ifdef DENDRO_PROFILE_UNZIP_DETAIL
profiler_t t_deriv_first;
profiler_t t_deriv_second;
long double t_ud_acc[UD_N];
static long double s_ud_before[UD_N];
static profiler_t* const* ud_timers() {
    namespace dt = dendro::timer;
    static profiler_t* const t[UD_N] = {
        &dt::t_ghost_pack,       &dt::t_ghost_wait,    &dt::t_ghost_unpack,
        &dt::t_ud_scatter_total, &dt::t_ud_nodal,      &dt::t_ud_same,
        &dt::t_ud_f2c,           &dt::t_ud_c2f_interp, &dt::t_ud_c2f_scat,
        &dt::t_ud_ring};
    return t;
}
void ud_snap_begin() {
    profiler_t* const* t = ud_timers();
    for (unsigned i = 0; i < UD_N; i++) s_ud_before[i] = t[i]->seconds;
}
void ud_snap_end() {
    profiler_t* const* t = ud_timers();
    for (unsigned i = 0; i < UD_N; i++)
        t_ud_acc[i] += t[i]->seconds - s_ud_before[i];
}
#endif

profiler_t t_rhs_a;
profiler_t t_rhs_b;
profiler_t t_rhs_gt;
profiler_t t_rhs_chi;
profiler_t t_rhs_At;
profiler_t t_rhs_K;
profiler_t t_rhs_Gt;
profiler_t t_rhs_B;

profiler_t t_bdyc;

profiler_t t_zip;
profiler_t t_rkStep;
profiler_t t_rkStage[6];

profiler_t t_isReMesh;
profiler_t t_gridTransfer;
profiler_t t_ioVtu;
profiler_t t_ioCheckPoint;

}  // namespace timer
}  // namespace bssn

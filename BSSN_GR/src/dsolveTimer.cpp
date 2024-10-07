
#include "dsolveTimer.h"

#include "dendroProfileParams.h"
#include "parameters.h"

namespace dsolve {

namespace timer {

profiler_t total_runtime;

profiler_t t_f2o;
profiler_t t_cons;
profiler_t t_bal;
profiler_t t_mesh;

profiler_t t_rkSolve;

profiler_t t_ghostEx_sync;

profiler_t t_unzip_sync;

profiler_t t_unzip_async;

profiler_t t_deriv;
profiler_t t_rhs;
profiler_t t_rhs_outer;

profiler_t t_bdyc;

profiler_t t_zip;
profiler_t t_rkStep;

profiler_t t_isReMesh;
profiler_t t_gridTransfer;
profiler_t t_ioVtu;
profiler_t t_ioCheckPoint;

void initFlops() {
    total_runtime.start();
    t_f2o.start();
    t_cons.start();
    t_bal.start();
    t_mesh.start();
    t_rkSolve.start();
    t_ghostEx_sync.start();
    t_unzip_sync.start();

    for (unsigned int i = 0; i < NUM_FACES; i++)
        dendro::timer::t_unzip_sync_face[i].start();

    dendro::timer::t_unzip_async_internal.start();
    dendro::timer::t_unzip_sync_edge.start();
    dendro::timer::t_unzip_sync_vtex.start();
    dendro::timer::t_unzip_p2c.start();
    dendro::timer::t_unzip_sync_nodalval.start();
    dendro::timer::t_unzip_sync_cpy.start();
    dendro::timer::t_unzip_sync_f_c1.start();
    dendro::timer::t_unzip_sync_f_c2.start();
    dendro::timer::t_unzip_sync_f_c3.start();

    t_unzip_async.start();
    dendro::timer::t_unzip_async_comm.start();

    dendro::timer::t_unzip_async_internal.start();
    dendro::timer::t_unzip_async_external.start();
    dendro::timer::t_unzip_async_comm.start();
    t_deriv.start();
    t_rhs.start();
    t_rhs_outer.start();

    dendro::timer::t_compression_extraction.start();
    dendro::timer::t_compression_compress.start();
    dendro::timer::t_compression_begin_comms.start();
    dendro::timer::t_compression_wait_comms.start();
    dendro::timer::t_compression_decompress.start();
    dendro::timer::t_compression_unextract.start();
    dendro::timer::t_compression_uzip_post.start();

    // t_rhs_a.start();
    // t_rhs_b.start();
    // t_rhs_gt.start();
    // t_rhs_chi.start();
    // t_rhs_At.start();
    // t_rhs_K.start();
    // t_rhs_Gt.start();
    // t_rhs_B.start();
    //
    t_bdyc.start();

    t_zip.start();
    t_rkStep.start();
    t_isReMesh.start();
    t_gridTransfer.start();
    t_ioVtu.start();
    t_ioCheckPoint.start();
}

void resetSnapshot() {
    total_runtime.snapreset();
    t_f2o.snapreset();
    t_cons.snapreset();
    t_bal.snapreset();
    t_mesh.snapreset();
    t_rkSolve.snapreset();
    t_ghostEx_sync.snapreset();
    t_unzip_sync.snapreset();

    for (unsigned int i = 0; i < NUM_FACES; i++)
        dendro::timer::t_unzip_sync_face[i].snapreset();

    dendro::timer::t_unzip_sync_internal.snapreset();
    dendro::timer::t_unzip_sync_edge.snapreset();
    dendro::timer::t_unzip_sync_vtex.snapreset();
    dendro::timer::t_unzip_p2c.snapreset();
    dendro::timer::t_unzip_sync_nodalval.snapreset();
    dendro::timer::t_unzip_sync_cpy.snapreset();

    dendro::timer::t_unzip_sync_f_c1.snapreset();
    dendro::timer::t_unzip_sync_f_c2.snapreset();
    dendro::timer::t_unzip_sync_f_c3.snapreset();

    t_unzip_async.snapreset();
    dendro::timer::t_unzip_async_internal.snapreset();
    dendro::timer::t_unzip_async_external.snapreset();
    dendro::timer::t_unzip_async_comm.snapreset();

    t_deriv.snapreset();
    t_rhs.snapreset();
    t_rhs_outer.snapreset();

    dendro::timer::t_compression_extraction.snapreset();
    dendro::timer::t_compression_compress.snapreset();
    dendro::timer::t_compression_begin_comms.snapreset();
    dendro::timer::t_compression_wait_comms.snapreset();
    dendro::timer::t_compression_decompress.snapreset();
    dendro::timer::t_compression_unextract.snapreset();
    dendro::timer::t_compression_uzip_post.snapreset();

    // t_rhs_a.snapreset();
    // t_rhs_b.snapreset();
    // t_rhs_gt.snapreset();
    // t_rhs_chi.snapreset();
    // t_rhs_At.snapreset();
    // t_rhs_K.snapreset();
    // t_rhs_Gt.snapreset();
    // t_rhs_B.snapreset();
    //
    t_bdyc.snapreset();

    t_zip.snapreset();
    t_rkStep.snapreset();
    t_isReMesh.snapreset();
    t_gridTransfer.snapreset();
    t_ioVtu.snapreset();
    t_ioCheckPoint.snapreset();
}

}  // namespace timer

}  // namespace dsolve

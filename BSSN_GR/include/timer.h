#pragma once

#include "bssnCtx.h"
#include "mesh.h"
#include "profiler.h"

namespace dsolve {

namespace timer {

extern profiler_t total_runtime;

extern profiler_t t_f2o;
extern profiler_t t_cons;
extern profiler_t t_bal;
extern profiler_t t_mesh;

extern profiler_t t_rkSolve;

extern profiler_t t_ghostEx_sync;

extern profiler_t t_unzip_sync;

extern profiler_t t_unzip_async;

extern profiler_t t_deriv;
extern profiler_t t_rhs;
extern profiler_t t_rhs_outer;

extern profiler_t t_bdyc;

extern profiler_t t_zip;
extern profiler_t t_rkStep;

extern profiler_t t_isReMesh;
extern profiler_t t_gridTransfer;
extern profiler_t t_ioVtu;
extern profiler_t t_ioCheckPoint;

/**@brief initialize all the flop counters. */
void initFlops();

/**@brief clears the snapshot counter for time profiler variables*/
void resetSnapshot();

/**@brief reduce min mean max.
 * @param [in] stat: local time
 * @param [out] stat_g 0-min, 1-mean 2-max
 * */
template <typename T>
void computeOverallStats(T *stat, T *stat_g, MPI_Comm comm) {
    int rank, npes;
    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &rank);

    par::Mpi_Reduce(stat, stat_g, 1, MPI_MIN, 0, comm);
    par::Mpi_Reduce(stat, stat_g + 1, 1, MPI_SUM, 0, comm);
    par::Mpi_Reduce(stat, stat_g + 2, 1, MPI_MAX, 0, comm);
    stat_g[1] /= (npes);
}

/** @breif : printout the profile parameters. */
void profileInfo(const char *filePrefix, const ot::Mesh *pMesh);

/** @breif : printout the profile parameters (intermediate profile information).
 */
template <typename CtxDerived>
void profileInfoIntermediate(
    const char *filePrefix, const ot::Mesh *pMesh,
    const unsigned int currentStep,
    const ts::Ctx<CtxDerived, DendroScalar, unsigned int> *ctxObj) {
    int activeRank, activeNpes, globalRank, globalNpes;

    MPI_Comm commActive;
    MPI_Comm commGlobal;

    if (pMesh->isActive()) {
        commActive = pMesh->getMPICommunicator();
        activeRank = pMesh->getMPIRank();
        activeNpes = pMesh->getMPICommSize();
    } else {
        return;
    }

    globalRank = pMesh->getMPIRankGlobal();
    globalNpes = pMesh->getMPICommSizeGlobal();
    commGlobal = pMesh->getMPIGlobalCommunicator();

    double t_stat;
    double t_stat_g[3];

    const char separator = ' ';
    const int nameWidth  = 30;
    const int numWidth   = 10;

    char fName[256];
    std::ofstream outfile;

    DendroIntL localSz, globalSz;

    DendroIntL ghostElements;
    DendroIntL localElements;

    DendroIntL ghostNodes;
    DendroIntL localNodes;

    DendroIntL totalSendNode;
    DendroIntL totalRecvNode;

    DendroIntL numCalls;

#ifdef SOLVER_PROFILE_HUMAN_READABLE
    if (!activeRank) {
        sprintf(fName, "%s_im.prof", filePrefix);
        outfile.open(fName, std::fstream::app);
        if (outfile.fail()) {
            std::cout << fName << " file open failed " << std::endl;
            return;
        }

        outfile << "active npes : " << activeNpes << std::endl;
        outfile << "global npes : " << globalNpes << std::endl;
        outfile << "current step : " << currentStep << std::endl;
        outfile << "partition tol : " << bssn::BSSN_LOAD_IMB_TOL << std::endl;
        outfile << "wavelet tol : " << bssn::BSSN_WAVELET_TOL << std::endl;
        outfile << "maxdepth : " << bssn::BSSN_MAXDEPTH << std::endl;
    }

    MPI_Comm comm     = commActive;
    unsigned int rank = activeRank;

    localSz           = pMesh->getNumLocalMeshElements();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << "Elements : " << globalSz << std::endl;

    localSz = pMesh->getNumLocalMeshNodes();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << "DOG(zip) : " << globalSz << std::endl;

    localSz = pMesh->getDegOfFreedomUnZip();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << "DOG(unzip) : " << globalSz << std::endl;

    ghostElements =
        pMesh->getNumPreGhostElements() + pMesh->getNumPostGhostElements();
    localElements = pMesh->getNumLocalMeshElements();

    ghostNodes    = pMesh->getNumPreMeshNodes() + pMesh->getNumPostMeshNodes();
    localNodes    = pMesh->getNumLocalMeshNodes();

    if (!rank)
        outfile << "========================= MESH "
                   "==========================================================="
                   "============ "
                << std::endl;

    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "step";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "min(#)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "mean(#)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "max(#)" << std::endl;

    t_stat = ghostElements;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "ghost Elements";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = localElements;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "local Elements";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = ghostNodes;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "ghost Nodes";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = localNodes;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "local Nodes";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = pMesh->getGhostExcgTotalSendNodeCount();
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "send Nodes";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = pMesh->getGhostExcgTotalRecvNodeCount();
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "recv Nodes";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    if (!rank)
        outfile << "========================= RUNTIME "
                   "==========================================================="
                   "======== "
                << std::endl;
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "step";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "min(s)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "mean(s)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "max(s)" << std::endl;

    /* t_stat=total_runtime.seconds;
    computeOverallStats(&t_stat,t_stat_g,comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<"+runtime(s)"; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


    t_stat=t_f2o.seconds;
    computeOverallStats(&t_stat,t_stat_g,comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<" ++f2o"; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


    t_stat=t_cons.seconds;
    computeOverallStats(&t_stat,t_stat_g,comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<" ++construction"; if(!rank)outfile << std::left
    << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


    t_stat=t_rkSolve.seconds;
    computeOverallStats(&t_stat,t_stat_g,comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<" ++rkSolve"; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;*/

    t_stat = t_bal.snap;
    // numCalls=t_bal.num_calls;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++2:1 balance";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_mesh.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++mesh";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rkStep.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++rkstep";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_ghostEx_sync.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++ghostExchge.";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_unzip_sync.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++unzip_sync";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_unzip_async.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++unzip_async";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

#ifdef ENABLE_DENDRO_PROFILE_COUNTERS

    t_stat = dendro::timer::t_unzip_async_comm.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_comm_wait (comm) ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_nodalval.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_nodalVal";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_f_c1.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --t_unzip_sync_f_c1";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_f_c2.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --t_unzip_sync_f_c2";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_f_c3.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --t_unzip_sync_f_c3";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_cpy.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --t_unzip_sync_cpy";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_internal.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_internal";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[0].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_left";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[1].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_right";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[2].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_down";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[3].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_up";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[4].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_back";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[5].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_front";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_edge.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_edge";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_vtex.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_vtex";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_p2c.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_p2c";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;
#endif

    /*
    #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
    t_stat=t_unzip_async_internal.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<"  --unzip_internal"; if(!rank)outfile <<
    std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

    t_stat=t_unzip_async_external.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<"  --unzip_external"; if(!rank)outfile <<
    std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


    t_stat=t_unzip_async_comm.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<"  --unzip_comm (comm) "; if(!rank)outfile <<
    std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;
    #endif
    */
    t_stat = t_isReMesh.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++isReMesh";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_gridTransfer.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++gridTransfer";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_deriv.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++deriv ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++compute_rhs ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    //    t_stat = t_rhs_a.snap;
    //    computeOverallStats(&t_stat, t_stat_g, comm);
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << "  --compute_rhs_a ";
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[0];
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[1];
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[2] << std::endl;
    //
    //    t_stat = t_rhs_b.snap;
    //    computeOverallStats(&t_stat, t_stat_g, comm);
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << "  --compute_rhs_b ";
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[0];
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[1];
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[2] << std::endl;
    //
    //    t_stat = t_rhs_gt.snap;
    //    computeOverallStats(&t_stat, t_stat_g, comm);
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << "  --compute_rhs_gt ";
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[0];
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[1];
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[2] << std::endl;
    //
    //    t_stat = t_rhs_chi.snap;
    //    computeOverallStats(&t_stat, t_stat_g, comm);
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << "  --compute_rhs_chi ";
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[0];
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[1];
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[2] << std::endl;
    //
    //    t_stat = t_rhs_At.snap;
    //    computeOverallStats(&t_stat, t_stat_g, comm);
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << "  --compute_rhs_At ";
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[0];
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[1];
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[2] << std::endl;
    //
    //    t_stat = t_rhs_K.snap;
    //    computeOverallStats(&t_stat, t_stat_g, comm);
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << "  --compute_rhs_K ";
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[0];
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[1];
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[2] << std::endl;
    //
    //    t_stat = t_rhs_Gt.snap;
    //    computeOverallStats(&t_stat, t_stat_g, comm);
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << "  --compute_rhs_Gt ";
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[0];
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[1];
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[2] << std::endl;
    //
    //    t_stat = t_rhs_B.snap;
    //    computeOverallStats(&t_stat, t_stat_g, comm);
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << "  --compute_rhs_B ";
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[0];
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[1];
    //    if (!rank)
    //        outfile << std::left << std::setw(nameWidth) <<
    //        std::setfill(separator)
    //                << t_stat_g[2] << std::endl;
    //
    t_stat = t_bdyc.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++boundary con ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_zip.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++zip";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_ioVtu.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++vtu";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_ioCheckPoint.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++checkpoint";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    if (!rank) outfile.close();
#else

    if (!activeRank) {
        sprintf(fName, "%s_im.prof", filePrefix);
        outfile.open(fName, std::fstream::app);
        if (outfile.fail()) {
            std::cout << fName << " file open failed " << std::endl;
            return;
        }

        // writes the header
        if (currentStep == 0)
            outfile
                << "step\t act_npes\t glb_npes\t part_tol\t wave_tol\t "
                   "maxdepth\t numOcts\t dof_zip\t dof_unzip\t"
                << "element_ghost_min\t element_ghost_mean\t "
                   "element_ghost_max\t"
                << "element_local_min\t element_local_mean\t "
                   "element_local_max\t"
                << "element_send_min\t element_send_mean\t element_send_max\t "
                << "element_recv_min\t element_recv_mean\t element_recv_max\t "
                << "element_send_min_comp\t element_send_mean_comp\t "
                   "element_send_max_comp\t "
                << "element_recv_min_comp\t element_recv_mean_comp\t "
                   "element_recv_max_comp\t "
                << "nodes_local_min\t nodes_local_mean\t nodes_local|max\t"
                << "send_nodes_min\t send_nodes_mean\t send_nodes_max\t"
                << "recv_nodes_min\t recv_nodes_mean\t recv_nodes_max\t"
                << "bal_min\t bal_mean\t bal_max\t"
                << "mesh_min\t mesh_mean\t mesh_max\t"
                << "rkstep_min\t rkstep_mean\t rkstep_max\t"
                << "ghostEx_min\t ghostEx_mean\t ghostEx_max\t"
                << "unzip_sync_min\t unzip_sync_mean\t unzip_sync_max\t"
                << "unzip_async_min\t unzip_async_mean\t unzip_async_max\t"
                << "unzip_async_wait_min\t unzip_async_wait_mean\t "
                   "unzip_async_wait_max\t"
                << "isRemesh_min\t isRemesh_mean\t isRemesh_max\t"
                << "GT_min\t GT_mean\t GT_max\t"
                << "zip_min\t zip_mean\t zip_max\t"
                << "deriv_min\t deriv_mean\t deriv_max\t"
                << "rhs_min\t rhs_mean\t rhs_max\t"
                << "bdyc_min\t bdyc_mean\t bdyc_max\t"
                << "rhs_outer_min\t rhs_outer_mean\t rhs_outer_max\t"
                << std::endl;
    }

    MPI_Comm comm     = commActive;
    unsigned int rank = activeRank;

    if (!rank) outfile << std::fixed << std::setprecision(15);

    if (!rank) outfile << currentStep << "\t ";
    if (!rank) outfile << activeNpes << "\t ";
    if (!rank) outfile << globalNpes << "\t ";
    if (!rank) outfile << bssn::BSSN_LOAD_IMB_TOL << "\t ";
    if (!rank) outfile << bssn::BSSN_WAVELET_TOL << "\t ";
    if (!rank) outfile << bssn::BSSN_MAXDEPTH << "\t ";

    localSz = pMesh->getNumLocalMeshElements();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << globalSz << "\t ";

    localSz = pMesh->getNumLocalMeshNodes();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << globalSz << "\t ";

    localSz = pMesh->getDegOfFreedomUnZip();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << globalSz << "\t ";

    ghostElements =
        pMesh->getNumPreGhostElements() + pMesh->getNumPostGhostElements();
    localElements = pMesh->getNumLocalMeshElements();

    t_stat        = ghostElements;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = localElements;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    // SEND AND RECEIVE DATA
    t_stat = (double)ctxObj->getTotalBytesSendSum();
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = (double)ctxObj->getTotalBytesRecvSum();
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = (double)ctxObj->getTotalBytesSendCompressSum();
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = (double)ctxObj->getTotalBytesRecvCompressSum();
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    ghostNodes = pMesh->getNumPreMeshNodes() + pMesh->getNumPostMeshNodes();
    localNodes = pMesh->getNumLocalMeshNodes();

    /*t_stat=ghostNodes;
    computeOverallStats(&t_stat,t_stat_g,comm);
    if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t
    ";*/

    t_stat     = localNodes;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = pMesh->getGhostExcgTotalSendNodeCount();
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = pMesh->getGhostExcgTotalRecvNodeCount();
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_bal.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_mesh.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_rkStep.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_ghostEx_sync.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_unzip_sync.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_unzip_async.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = dendro::timer::t_unzip_async_comm.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_isReMesh.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_gridTransfer.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_zip.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_deriv.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_rhs.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_bdyc.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_rhs_outer.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    if (!rank) outfile << std::endl;
    if (!rank) outfile.close();
#endif
}

}  // namespace timer

}  // namespace dsolve

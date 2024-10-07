/**
 * @file gr_scaling.cpp
 * @brief BSSN_GR scaling tests.
 * @version 0.1
 * @date 2021-07-25
 *
 * @copyright Copyright (c) 2021
 *
 */

#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>

#include "TreeNode.h"
#include "bssnCtx.h"
#include "dsolveTimer.h"
#include "gr.h"
#include "grUtils.h"
#include "mathUtils.h"
#include "mesh.h"
#include "meshUtils.h"
#include "mpi.h"
#include "octUtils.h"
#include "parameters.h"

int bssn_driver(MPI_Comm comm, unsigned int num_step, unsigned int warm_up,
                std::ostream& outfile, unsigned int ts_mode) {
    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double, double, double, double*)> f_init =
        [](double x, double y, double z, double* var) {
            bssn::punctureData(x, y, z, var);
        };
    const unsigned int interpVars = bssn::BSSN_NUM_VARS;

    unsigned int varIndex[interpVars];
    for (unsigned int i = 0; i < bssn::BSSN_NUM_VARS; i++) varIndex[i] = i;

    if (false && bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY) {
        if (!rank)
            std::cout << YLW << "Using block adaptive mesh. AMR disabled "
                      << NRM << std::endl;
        const Point pt_min(bssn::BSSN_BLK_MIN_X, bssn::BSSN_BLK_MIN_Y,
                           bssn::BSSN_BLK_MIN_Z);
        const Point pt_max(bssn::BSSN_BLK_MAX_X, bssn::BSSN_BLK_MAX_Y,
                           bssn::BSSN_BLK_MAX_Z);
        bssn::blockAdaptiveOctree(tmpNodes, pt_min, pt_max, m_uiMaxDepth - 2,
                                  m_uiMaxDepth, comm);
    } else {
        if (!rank)
            std::cout << YLW << "Using function2Octree. AMR enabled " << NRM
                      << std::endl;
        const unsigned int f2olmin =
            std::min(bssn::BSSN_BH1_MAX_LEV, bssn::BSSN_BH2_MAX_LEV);
        if (f2olmin < MAXDEAPTH_LEVEL_DIFF + 2) {
            if (!rank)
                std::cout << "BH min level should be larger than "
                          << (MAXDEAPTH_LEVEL_DIFF + 2) << std::endl;

            MPI_Abort(comm, 0);
        }
        function2Octree(f_init, bssn::BSSN_NUM_VARS, varIndex, interpVars,
                        tmpNodes, (f2olmin - MAXDEAPTH_LEVEL_DIFF - 2),
                        bssn::BSSN_WAVELET_TOL, bssn::BSSN_ELE_ORDER, comm);
    }

    // std::vector<ot::TreeNode> f2Octants(tmpNodes);
    ot::Mesh* mesh =
        ot::createMesh(tmpNodes.data(), tmpNodes.size(), bssn::BSSN_ELE_ORDER,
                       comm, 1, ot::SM_TYPE::FDM, bssn::BSSN_DENDRO_GRAIN_SZ,
                       bssn::BSSN_LOAD_IMB_TOL, bssn::BSSN_SPLIT_FIX);
    mesh->setDomainBounds(Point(bssn::BSSN_GRID_MIN_X, bssn::BSSN_GRID_MIN_Y,
                                bssn::BSSN_GRID_MIN_Z),
                          Point(bssn::BSSN_GRID_MAX_X, bssn::BSSN_GRID_MAX_Y,
                                bssn::BSSN_GRID_MAX_Z));
    unsigned int lmin, lmax;
    mesh->computeMinMaxLevel(lmin, lmax);
    tmpNodes.clear();

    if (!rank) {
        std::cout << "================= Grid Info (Before init grid "
                     "converge):==============================================="
                     "========"
                  << std::endl;
        std::cout << "lmin: " << lmin << " lmax:" << lmax << std::endl;
        std::cout << "dx: "
                  << ((bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
                      ((1u << (m_uiMaxDepth - lmax)) /
                       ((double)bssn::BSSN_ELE_ORDER)) /
                      ((double)(1u << (m_uiMaxDepth))))
                  << std::endl;
        std::cout << "dt: "
                  << bssn::BSSN_CFL_FACTOR *
                         ((bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
                          ((1u << (m_uiMaxDepth - lmax)) /
                           ((double)bssn::BSSN_ELE_ORDER)) /
                          ((double)(1u << (m_uiMaxDepth))))
                  << std::endl;
        std::cout << "ts mode: " << ts_mode << std::endl;
        std::cout << "========================================================="
                     "======================================================"
                  << std::endl;
    }

    // calculate the minimum dx
    bssn::BSSN_CURRENT_MIN_DX =
        ((bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
         ((1u << (m_uiMaxDepth - lmax)) / ((double)bssn::BSSN_ELE_ORDER)) /
         ((double)(1u << (m_uiMaxDepth))));

    bssn::BSSN_RK45_TIME_STEP_SIZE =
        bssn::BSSN_CFL_FACTOR *
        ((bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
         ((1u << (m_uiMaxDepth - lmax)) / ((double)bssn::BSSN_ELE_ORDER)) /
         ((double)(1u << (m_uiMaxDepth))));
    tmpNodes.clear();

    // enable block adaptivity disable remesh.
    bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY = 1;
    ot::Mesh* pMesh                    = mesh;

    if (bssn::BSSN_RESTORE_SOLVER == 0)
        pMesh = bssn::weakScalingReMesh(mesh, npes);

    if (ts_mode == 0) {
    } else if (ts_mode == 1) {
        bssn::BSSNCtx* bssnCtx = new bssn::BSSNCtx(pMesh);
        ts::ETS<DendroScalar, bssn::BSSNCtx>* ets =
            new ts::ETS<DendroScalar, bssn::BSSNCtx>(bssnCtx);
        ets->set_evolve_vars(bssnCtx->get_evolution_vars());

        if ((RKType)bssn::BSSN_RK_TYPE == RKType::RK3)
            ets->set_ets_coefficients(ts::ETSType::RK3);
        else if ((RKType)bssn::BSSN_RK_TYPE == RKType::RK4)
            ets->set_ets_coefficients(ts::ETSType::RK4);
        else if ((RKType)bssn::BSSN_RK_TYPE == RKType::RK45)
            ets->set_ets_coefficients(ts::ETSType::RK5);

        ets->init();

        for (unsigned int ii = 0; ii < (warm_up + num_step); ii++) {
            const DendroIntL step   = ets->curr_step();
            const DendroScalar time = ets->curr_time();
            pMesh                   = bssnCtx->get_mesh();

#if defined __PROFILE_ETS__ && __PROFILE_CTX__
            if (step == warm_up) {
                ets->reset_pt();
                bssnCtx->reset_pt();
            }
#endif

            if (step == 0) bssnCtx->write_checkpt();

            bssn::BSSN_CURRENT_RK_COORD_TIME = step;
            bssn::BSSN_CURRENT_RK_STEP       = time;

            const bool isActive              = ets->is_active();
            const unsigned int rank_global   = ets->get_global_rank();

            if (!rank_global)
                std::cout << "[ETS] : Executing step :  " << ets->curr_step()
                          << "\tcurrent time :" << ets->curr_time()
                          << "\t dt:" << ets->ts_size() << "\t" << std::endl;

            bssnCtx->resetCountersForEvolve();

            dsolve::timer::t_rkStep.start();
            ets->evolve();
            dsolve::timer::t_rkStep.stop();

            bssnCtx->resetForNextStep();

            MPI_Barrier(comm);
        }

        char fname[256];
        sprintf(fname, "bssnCtx_WS_%d_", npes);

        dsolve::timer::profileInfoIntermediate(fname, bssnCtx->get_mesh(),
                                               ets->curr_step(), bssnCtx);

        ot::Mesh* tmp_mesh = bssnCtx->get_mesh();
        delete bssnCtx;
        delete tmp_mesh;
        delete ets;
    } else {
        if (!rank) RAISE_ERROR("invalid ts mode : " << ts_mode << "specifed");

        MPI_Abort(comm, 0);
    }

    return 0;
}

int main(int argc, char** argv) {
    // 0- NUTS 1-UTS
    unsigned int ts_mode = 0;

    if (argc < 2) {
        std::cout << "Usage: " << argv[0]
                  << "paramFile TSMode(0){0-Spatially Adaptive Time "
                     "Stepping(SATS, "
                  << GRN << "default" << NRM
                  << ") , 1- Uniform Time Stepping.  }" << std::endl;
        return 0;
    }

    if (argc > 2) ts_mode = std::atoi(argv[2]);

    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    // Print out CMAKE options
    if (!rank) {
#ifdef BSSN_COMPUTE_CONSTRAINTS
        std::cout << GRN << "  Compiled with BSSN_COMPUTE_CONSTRAINTS" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without BSSN_COMPUTE_CONSTRAINTS" << NRM
                  << std::endl;
#endif
#ifdef BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT
        std::cout << GRN << "  Compiled with BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT"
                  << NRM << std::endl;
#else
        std::cout << RED
                  << "  Compiled without BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT"
                  << NRM << std::endl;
#endif
#ifdef BSSN_ENABLE_VTU_OUTPUT
        std::cout << GRN << "  Compiled with BSSN_ENABLE_VTU_OUTPUT" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without BSSN_ENABLE_VTU_OUTPUT" << NRM
                  << std::endl;
#endif
#ifdef BSSN_ETA_FUNCTION
        std::cout << GRN << "  Compiled with  BSSN_ETA_FUNCTION" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without  BSSN_ETA_FUNCTION" << NRM
                  << std::endl;
#endif
#ifdef BSSN_EXTRACT_BH_LOCATIONS
        std::cout << GRN << "  Compiled with  BSSN_EXTRACT_BH_LOCATIONS" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without  BSSN_EXTRACT_BH_LOCATIONS"
                  << NRM << std::endl;
#endif
#ifdef BSSN_EXTRACT_GRAVITATIONAL_WAVES
        std::cout << GRN << "  Compiled with  BSSN_EXTRACT_GRAVITATIONAL_WAVES"
                  << NRM << std::endl;
#else
        std::cout << RED
                  << "  Compiled without  BSSN_EXTRACT_GRAVITATIONAL_WAVES"
                  << NRM << std::endl;
#endif
#ifdef BSSN_EXTRACT_GRAVITATIONAL_WAVES
        std::cout << GRN << "  Compiled with  BSSN_EXTRACT_GRAVITATIONAL_WAVES"
                  << NRM << std::endl;
#else
        std::cout << RED
                  << "  Compiled without  BSSN_EXTRACT_GRAVITATIONAL_WAVES"
                  << NRM << std::endl;
#endif
#ifdef BSSN_GAUGE_ROCHESTER
        std::cout << GRN << "  Compiled with  BSSN_GAUGE_ROCHESTER" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without  BSSN_GAUGE_ROCHESTER" << NRM
                  << std::endl;
#endif
#ifdef BSSN_KERR_SCHILD_TEST
        std::cout << GRN << "  Compiled with  BSSN_KERR_SCHILD_TEST" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without  BSSN_KERR_SCHILD_TEST" << NRM
                  << std::endl;
#endif

#ifdef BSSN_REFINE_BASE_EH
        std::cout << GRN << "  Compiled with  BSSN_REFINE_BASE_EH" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without  BSSN_REFINE_BASE_EH" << NRM
                  << std::endl;
#endif

#ifdef USE_FD_INTERP_FOR_UNZIP
        std::cout << GRN << "  Compiled with  USE_FD_INTERP_FOR_UNZIP" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without  USE_FD_INTERP_FOR_UNZIP" << NRM
                  << std::endl;
#endif
    }

    // 1 . read the parameter file.
    if (!rank) std::cout << " reading parameter file :" << argv[1] << std::endl;
    bssn::readParamFile(argv[1], comm);

    int root = std::min(1, npes - 1);
    bssn::dumpParamFile(std::cout, root, comm);

    _InitializeHcurve(bssn::BSSN_DIM);
    m_uiMaxDepth = bssn::BSSN_MAXDEPTH;

    if (bssn::BSSN_NUM_VARS % bssn::BSSN_ASYNC_COMM_K != 0) {
        if (!rank)
            std::cout << "[overlap communication error]: total BSSN_NUM_VARS: "
                      << bssn::BSSN_NUM_VARS
                      << " is not divisable by BSSN_ASYNC_COMM_K: "
                      << bssn::BSSN_ASYNC_COMM_K << std::endl;
        MPI_Abort(comm, 0);
    }

    if (bssn::BSSN_GW_EXTRACT_FREQ > bssn::BSSN_IO_OUTPUT_FREQ) {
        if (!rank)
            std::cout << " BSSN_GW_EXTRACT_FREQ  should be less "
                         "BSSN_IO_OUTPUT_FREQ "
                      << std::endl;
        MPI_Abort(comm, 0);
    }

    // 2. generate the initial grid.
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double, double, double, double*)> f_init =
        [](double x, double y, double z, double* var) {
            bssn::punctureData(x, y, z, var);
        };
    std::function<double(double, double, double)> f_init_alpha =
        [](double x, double y, double z) {
            double var[24];
            bssn::punctureData(x, y, z, var);
            return var[0];
        };
    // std::function<void(double,double,double,double*)> f_init=[](double
    // x,double y,double z,double*var){bssn::KerrSchildData(x,y,z,var);};

    const unsigned int interpVars = bssn::BSSN_NUM_VARS;
    unsigned int varIndex[interpVars];
    for (unsigned int i = 0; i < bssn::BSSN_NUM_VARS; i++) varIndex[i] = i;

    /*varIndex[0]=bssn::VAR::U_ALPHA;
    varIndex[1]=bssn::VAR::U_CHI;*/
    DendroIntL localSz, globalSz;
    double t_stat;
    double t_stat_g[3];

    const unsigned int NUM_WARM_UP = 2;
    const unsigned int NUM_STEPS   = 4;

    std::ofstream outfile;
    char fname[256];
    sprintf(fname, "bssnCtx_WS_%d.txt", npes);

    if (!rank) {
        outfile.open(fname, std::ios_base::app);
        time_t now = time(0);
        // convert now to string form
        char* dt   = ctime(&now);
        outfile << "======================================================="
                   "====="
                << std::endl;
        outfile << "Current time : " << dt << " --- " << std::endl;
        outfile << "======================================================="
                   "====="
                << std::endl;
    }

    bssn_driver(comm, NUM_STEPS, NUM_WARM_UP, outfile, 1);

    if (!rank) outfile.close();

#ifdef RUN_WEAK_SCALING

    if (!rank)
        std::cout << "========================================================="
                     "============="
                  << std::endl;
    if (!rank) std::cout << "     Weak Scaling Run Begin.   " << std::endl;
    if (!rank)
        std::cout << "========================================================="
                     "============="
                  << std::endl;

    int proc_group = 0;
    int min_np     = 2;
    for (int i = npes; rank < i && i >= min_np; i = i >> 1) proc_group++;
    MPI_Comm comm_ws;

    MPI_Comm_split(comm, proc_group, rank, &comm_ws);

    MPI_Comm_rank(comm_ws, &rank);
    MPI_Comm_size(comm_ws, &npes);

    if (!rank) outfile.open(fname, std::ios_base::app);
    MPI_Barrier(comm_ws);

    bssn_driver(comm_ws, NUM_STEPS, NUM_WARM_UP, outfile, 1);

    MPI_Barrier(comm_ws);
    if (!rank) outfile.close();

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);
    MPI_Barrier(comm);

    if (!rank)
        std::cout << "========================================================="
                     "============="
                  << std::endl;
    if (!rank) std::cout << "     Weak Scaling Run Complete.   " << std::endl;
    if (!rank)
        std::cout << "========================================================="
                     "============="
                  << std::endl;

#endif

    MPI_Finalize();
    return 0;
}

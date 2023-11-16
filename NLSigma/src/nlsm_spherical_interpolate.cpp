/**
 * @file nlsm_nuts.cpp
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief : NLSM to test with NUTS.
 * @version 0.1
 * @date 2020-04-03
 * @copyright Copyright (c) 2020
 *
 */

#include <cstdint>
#define NLSM_FORCE_SAVE_BEFORE_AND_AFTER

#include <climits>
#include <iostream>
#include <vector>

#include "TreeNode.h"
#include "assert.h"
#include "enuts.h"
#include "ets.h"
#include "mathUtils.h"
#include "mesh.h"
#include "mpi.h"
#include "nlsm.h"
#include "nlsmCtx.h"
#include "nlsmUtils.h"
#include "octUtils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

std::vector<double> sphereical_points(double xcenter, double ycenter,
                                      double zcenter, double r,
                                      uint32_t npoints) {
    // the output vector is a straight list of points of 3's
    std::vector<double> out_vector;
    // out_vector.resize(3 * npoints);

    // sphere equation is (x - xcenter)^2 + (y - ycenter)^2 + (z - zcenter)^2 =
    // r^2
    //
    // this is the fibonacci sphere algorithm that generates a uniform shell

    // golden angle in radians
    double phi = M_PI * (std::sqrt(5.0) - 1.0);

    double x, y, z, radius, theta;

    // coords for spherical location
    double r_t, th, ph;

    for (uint32_t ii = 0; ii < npoints; ii++) {
        // this is y from 1 to -1
        y = 1.0 - ((double)ii / ((double)npoints - 1)) * 2.0;
        // calculate y's radius
        radius = std::sqrt(1.0 - y * y);
        // adjust the golden angle
        theta = phi * (double)ii;

        // then calculate x and z
        x = std::cos(theta) * radius;
        z = std::sin(theta) * radius;

        // these are now a unit sphere centered at 0, so we'll want to shift
        // them out first adjust the radius

        r_t = std::sqrt(x * x + y * y + z * z);
        th = std::acos(z / r_t);
        ph = std::atan2(y, x);

        // now we can adjust our values to go back to x, y, and z, and move by
        // center
        x = r * std::sin(th) * std::cos(ph) + xcenter;
        y = r * std::sin(th) * std::sin(ph) + ycenter;
        z = r * std::cos(th) + zcenter;

        out_vector.push_back(x);
        out_vector.push_back(y);
        out_vector.push_back(z);
    }

    return out_vector;
}

// std::vector<double> spherical_via_sunflower(double xcenter, double ycenter,
//                                             double zcenter, double r,
//                                             uint32_t npoints) {
//     std::vector<double> out_vector;
//
//     return out_vector;
// }

int main(int argc, char** argv) {
    if (argc < 2)
        std::cout << "Usage: " << argv[0] << " paramFile" << std::endl;

    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    if (!rank) {
#ifdef NLSM_NONLINEAR
        std::cout << GRN << "Compiled with NLSM_NONLINEAR" << NRM << std::endl;
#else
        std::cout << RED << "Compiled without NLSM_NONLINEAR" << NRM
                  << std::endl;
#endif

#ifdef NLSM_COMPARE_WITH_ANALYTICAL_SOL
        std::cout << GRN << "Compiled with NLSM_COMPARE_WITH_ANALYTICAL_SOL"
                  << NRM << std::endl;
#else
        std::cout << RED << "Compiled without NLSM_COMPARE_WITH_ANALYTICAL_SOL"
                  << NRM << std::endl;
#endif

#ifdef NLSM_USE_4TH_ORDER_DERIVS
        std::cout << GRN << "Compiled with NLSM_USE_4TH_ORDER_DERIVS" << NRM
                  << std::endl;
#else
        std::cout << RED << "Compiled without NLSM_USE_4TH_ORDER_DERIVS" << NRM
                  << std::endl;
#endif

#ifdef NLSM_USE_6TH_ORDER_DERIVS
        std::cout << GRN << "Compiled with NLSM_USE_6TH_ORDER_DERIVS" << NRM
                  << std::endl;
#else
        std::cout << RED << "Compiled without NLSM_USE_6TH_ORDER_DERIVS" << NRM
                  << std::endl;
#endif

#ifdef NLSM_USE_8TH_ORDER_DERIVS
        std::cout << GRN << "Compiled with NLSM_USE_8TH_ORDER_DERIVS" << NRM
                  << std::endl;
#else
        std::cout << RED << "Compiled without NLSM_USE_8TH_ORDER_DERIVS" << NRM
                  << std::endl;
#endif
    }

    nlsm::timer::initFlops();

    nlsm::timer::total_runtime.start();

    /**
     * ==1== READ IN THE PARAMETER FILE
     *
     */
    if (!rank) std::cout << " reading parameter file :" << argv[1] << std::endl;
    nlsm::readParamFile(argv[1], comm);
    nlsm::dumpParamFile(std::cout, 1, comm);

    // initialize HCurve data, used for interpolation
    _InitializeHcurve(nlsm::NLSM_DIM);
    m_uiMaxDepth = nlsm::NLSM_MAXDEPTH;

    if (nlsm::NLSM_NUM_VARS % nlsm::NLSM_ASYNC_COMM_K != 0) {
        if (!rank)
            std::cout << "[overlap communication error]: total NLSM_NUM_VARS: "
                      << nlsm::NLSM_NUM_VARS
                      << " is not divisable by NLSM_ASYNC_COMM_K: "
                      << nlsm::NLSM_ASYNC_COMM_K << std::endl;
        exit(0);
    }

    /**
     * ==2== GENERATE INITIAL GRID
     *
     * Start by defining the initial data function, and then build up the grid
     * based on configurationns
     */
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double, double, double, double*)> f_init =
        [](double x, double y, double z, double* var) {
            nlsm::initData(x, y, z, var);
        };
    std::function<void(double, double, double, double, double*)> u_x_t =
        [](double x, double y, double z, double t, double* var) {
            nlsm::analyticalSol(x, y, z, t, var);
        };
    // std::function<void(double,double,double,double*)> f_init=[](double
    // x,double y,double z,double*var){nlsm::KerrSchildData(x,y,z,var);};

    const unsigned int interpVars = nlsm::NLSM_NUM_VARS;
    unsigned int varIndex[interpVars];
    for (unsigned int i = 0; i < nlsm::NLSM_NUM_VARS; i++) varIndex[i] = i;

    DendroIntL localSz, globalSz;
    double t_stat;
    double t_stat_g[3];

    // generate the intiial grid depending on if it's block adaptivity or
    // function2Octree
    if (nlsm::NLSM_ENABLE_BLOCK_ADAPTIVITY) {
        // Block Adaptivity doesn't do any adaptive refinement, it's useful for
        // setting a region of "importance" that it will refine over. if the
        // entire compute region is set as important, it will build a uniform
        // grid
        if (!rank)
            std::cout << YLW << "Using block adaptive mesh. AMR disabled "
                      << NRM << std::endl;
        const Point pt_min(nlsm::NLSM_BLK_MIN_X, nlsm::NLSM_BLK_MIN_Y,
                           nlsm::NLSM_BLK_MIN_Z);
        const Point pt_max(nlsm::NLSM_BLK_MAX_X, nlsm::NLSM_BLK_MAX_Y,
                           nlsm::NLSM_BLK_MAX_Z);

        nlsm::blockAdaptiveOctree(
            tmpNodes, pt_min, pt_max,
            m_uiMaxDepth - (binOp::fastLog2(nlsm::NLSM_ELE_ORDER)),
            m_uiMaxDepth, comm);
    } else {
        // block adaptive off means function2Octree. This will build a tree
        // based on the function above initData. initData on its own isn't
        // that useful in this function, but it is a good starting point
        if (!rank)
            std::cout << YLW << "Using function2Octree. AMR enabled " << NRM
                      << std::endl;
        function2Octree(f_init, nlsm::NLSM_NUM_VARS,
                        nlsm::NLSM_REFINE_VARIABLE_INDICES,
                        nlsm::NLSM_NUM_REFINE_VARS, tmpNodes, m_uiMaxDepth,
                        nlsm::NLSM_WAVELET_TOL, nlsm::NLSM_ELE_ORDER, comm);
    }
    // at this point we will have a built up tree

    unsigned int lmin = 1, lmax = 5;

    // Now build up a mesh object from the nodes themselves. The mesh contains
    // a lot more information and functionality!
    ot::Mesh* mesh =
        ot::createMesh(tmpNodes.data(), tmpNodes.size(), nlsm::NLSM_ELE_ORDER,
                       comm, 1, ot::SM_TYPE::FDM, nlsm::NLSM_DENDRO_GRAIN_SZ,
                       nlsm::NLSM_LOAD_IMB_TOL, nlsm::NLSM_SPLIT_FIX);
    mesh->setDomainBounds(Point(nlsm::NLSM_GRID_MIN_X, nlsm::NLSM_GRID_MIN_Y,
                                nlsm::NLSM_GRID_MIN_Z),
                          Point(nlsm::NLSM_GRID_MAX_X, nlsm::NLSM_GRID_MAX_Y,
                                nlsm::NLSM_GRID_MAX_Z));

    // This performs a refinement for us, in the BSSN solver this has been moved
    // to ctx object construction
    bool is_mindepth_refine_g = false;
    do {
        if (!rank)
            std::cout << "enforce min depth refinement currently only works "
                         "for block AMR for NLSM"
                      << std::endl;

        bool is_mindepth_refine = false;
        std::vector<unsigned int> refine_flag;
        refine_flag.reserve(mesh->getNumLocalMeshElements());
        const ot::TreeNode* pNodes = mesh->getAllElements().data();
        for (unsigned int ele = mesh->getElementLocalBegin();
             ele < mesh->getElementLocalEnd(); ele++) {
            if (pNodes[ele].getLevel() < nlsm::NLSM_MINDEPTH) {
                refine_flag.push_back(OCT_SPLIT);
                is_mindepth_refine = true;
            } else {
                refine_flag.push_back(OCT_NO_CHANGE);
            }
        }

        MPI_Allreduce(&is_mindepth_refine, &is_mindepth_refine_g, 1, MPI_C_BOOL,
                      MPI_LOR, comm);

        if (is_mindepth_refine_g) {
            mesh->setMeshRefinementFlags(refine_flag);
            ot::Mesh* newMesh = mesh->ReMesh();

            DendroIntL localSz = mesh->getNumLocalMeshElements();
            DendroIntL gSz_new, gSz_old;

            par::Mpi_Reduce(&localSz, &gSz_old, 1, MPI_SUM, 0, comm);
            localSz = newMesh->getNumLocalMeshElements();
            par::Mpi_Reduce(&localSz, &gSz_new, 1, MPI_SUM, 0, comm);

            if (!rank)
                std::cout << "old mesh size: " << gSz_old
                          << " new mesh size: " << gSz_new << std::endl;

            std::swap(newMesh, mesh);
            delete newMesh;
        }

    } while (is_mindepth_refine_g);

    mesh->computeMinMaxLevel(lmin, lmax);
    nlsm::NLSM_RK45_TIME_STEP_SIZE =
        nlsm::NLSM_CFL_FACTOR *
        ((nlsm::NLSM_COMPD_MAX[0] - nlsm::NLSM_COMPD_MIN[0]) *
         ((1u << (m_uiMaxDepth - lmax)) / ((double)nlsm::NLSM_ELE_ORDER)) /
         ((double)(1u << (m_uiMaxDepth))));
    par::Mpi_Bcast(&nlsm::NLSM_RK45_TIME_STEP_SIZE, 1, 0, comm);

    // and now we can print some information about the grid after it was refined
    // slightly
    DendroIntL lblocks = mesh->getLocalBlockList().size();
    DendroIntL gblocks = 0;
    par::Mpi_Reduce(&lblocks, &gblocks, 1, MPI_SUM, 0, comm);
    if (!rank)
        std::cout << " number of blocks for coarset block level : "
                  << (m_uiMaxDepth - MAXDEAPTH_LEVEL_DIFF - 1)
                  << " # blocks: " << gblocks << std::endl;

    if (!rank) std::cout << " lmin: " << lmin << " lmax: " << lmax << std::endl;

    if (!rank) std::cout << "ts_mode: " << ts_mode << std::endl;

    const ts::ETSType tsType = ts::ETSType::RK4;

    /**
     * ==4== CREATE CTX OBJECT
     *
     * The CTX object is what will store the mesh and provide a lot more power
     * regarding the individual variables. It has all of the data vectors
     */
    // UTS
    nlsm::NLSMCtx* appCtx = new nlsm::NLSMCtx(mesh);
    ts::ETS<DendroScalar, nlsm::NLSMCtx>* ets =
        new ts::ETS<DendroScalar, nlsm::NLSMCtx>(appCtx);
    ets->set_evolve_vars(appCtx->get_evolution_vars());
    ets->set_ets_coefficients(tsType);

    for (ets->init(); ets->curr_time() < nlsm::NLSM_RK45_TIME_END;
         ets->evolve()) {
        const DendroIntL step = ets->curr_step();
        const DendroScalar time = ets->curr_time();

        const bool isActive = ets->is_active();
        const unsigned int rank_global = ets->get_global_rank();

        if ((step % nlsm::NLSM_REMESH_TEST_FREQ) == 0) {
            bool isRemesh = appCtx->is_remesh();
            if (isRemesh) {
                if (!rank_global)
                    std::cout << "[ETS] : Remesh is triggered.  \n";

                appCtx->remesh_and_gridtransfer(nlsm::NLSM_DENDRO_GRAIN_SZ,
                                                nlsm::NLSM_LOAD_IMB_TOL,
                                                nlsm::NLSM_SPLIT_FIX);
                appCtx->terminal_output();

                ets->sync_with_mesh();
            }

            if ((step % nlsm::NLSM_IO_OUTPUT_FREQ) == 0) appCtx->write_vtu();

            if ((step % nlsm::NLSM_CHECKPT_FREQ) == 0) appCtx->write_checkpt();
        }
    }

#ifdef NLSM_FORCE_SAVE_BEFORE_AND_AFTER
    // so we can view a before and after...
    unsigned int prev_timestep =
        appCtx->force_update_current_step(UINT_MAX - 1);

    // OPTIONAL: write to a vtu the
    appCtx->write_vtu();
    // then go to UINT_MAX for our after
    unsigned int garbo_value = appCtx->force_update_current_step(UINT_MAX);
#endif

    // TODO: THIS IS WHERE WE DO THE MAGIC!
    //
    //
    //
    //
    //
    //
    //
    //

#ifdef NLSM_FORCE_SAVE_BEFORE_AND_AFTER
    // OPTIONAL: write to a vtu the "fixed" output
    appCtx->write_vtu();
    // return the appCTX object back to our previous timestep value for the vtu
    garbo_value = appCtx->force_update_current_step(prev_timestep);
#endif

    // IMPORTANT: we have to write the checkpoint **back** out to save all
    // changes
    appCtx->write_checkpt();

    delete appCtx->get_mesh();
    delete appCtx;
    delete ets;

    MPI_Finalize();

    return 0;
}

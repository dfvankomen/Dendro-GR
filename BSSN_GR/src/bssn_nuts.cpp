/**
 * @file bssn_nuts.cpp
 * @brief BSSN driver for the ENUTS (explicit non-uniform / local time stepper),
 *  ported from em4_nuts.cpp (which was ported from nlsm_nuts.cpp). BSSN/BBH is
 *  the real production target for LTS, so this is the vehicle for measuring
 *  whether LTS preserves the RK4 convergence order on the actual NONLINEAR
 *  system.
 *
 *  ts_mode 0 = ENUTS (LTS, sub-cycling)      [default]
 *  ts_mode 1 = ETS   (GTS, global uniform RK) -- the same integrator bssnSolver
 *              runs, so it is the trusted same-scheme reference / RK4 control.
 *
 *  Accuracy protocol (see the Repartitioning "Next Session Handoff"):
 *   BBH has NO analytical solution, so accuracy is measured by SELF-CONVERGENCE:
 *     p = log2( ||u(dt) - u(dt/2)|| / ||u(dt/2) - u(dt/4)|| )
 *   on a FIXED mesh (remesh off), SFC (BSSN_PARTITIONING_METHOD=0),
 *   DENDRO_NUTS_NO_WPART=1 (keeps ENUTS on the same grid as GTS), aligned final
 *   times. GTS self-convergence must come out = 4.00 (the RK4 control that
 *   validates the pipeline). See Logs/2026-07-17 Session 3 for the EM4 result.
 *
 *  Do NOT use the EM4 SOLVER_ETA_R0 trick here: BSSN's eta is the physical
 *  Gamma-driver gauge damping (feeds the RHS). maxdepth-11 punctures span many
 *  levels, so the natural (eta-limited) offset still sub-cycles heavily -- just
 *  run naturally and confirm `est speedup > 1`.
 */

#include <iomanip>
#include <iostream>
#include <vector>

#include "TreeNode.h"
#include "bssnCtx.h"
#include "gr.h"       // pulls in enuts.h (ts::ExplicitNUTS) + ets.h (ts::ETS)
#include "grUtils.h"
#include "logger.h"
#include "mesh.h"
#include "meshUtils.h"
#include "mpi.h"
#include "octUtils.h"
#include "parameters.h"

// DENDRO_EVAR_DUMP=<prefix>: write this rank's local evolution vector to
// <prefix>_rank<globalRank>.bin. Two separate runs (ENUTS ts_mode 0 w/
// DENDRO_NUTS_NO_WPART=1, GTS ts_mode 1) on the SAME fixed mesh (remesh off,
// same octree + npes -> deterministic SFC partition) produce byte-aligned files
// per rank -> diff offline = pointwise LTS-vs-GTS error. Process-isolated, so it
// sidesteps the MPI comm-lifecycle issues of two meshes in one process.
static void dump_evar_bin(const ot::Mesh* mesh, bssn::DVec& evar) {
    const char* pfx = std::getenv("DENDRO_EVAR_DUMP");
    if (!pfx || !mesh->isActive()) return;
    const unsigned int nlb  = mesh->getNodeLocalBegin();
    const unsigned int nloc = mesh->getNumLocalMeshNodes();
    const unsigned int dof  = evar.get_dof();
    DendroScalar* p[bssn::BSSN_NUM_VARS];
    evar.to_2d(p);
    char fn[512];
    snprintf(fn, sizeof(fn), "%s_rank%d.bin", pfx, mesh->getMPIRankGlobal());
    FILE* f = fopen(fn, "wb");
    if (!f) return;
    fwrite(&nloc, sizeof(unsigned int), 1, f);
    fwrite(&dof, sizeof(unsigned int), 1, f);
    for (unsigned int v = 0; v < dof; v++)
        fwrite(p[v] + nlb, sizeof(DendroScalar), nloc, f);
    fclose(f);
}

// After an SFC remesh the base ts::Ctx::remesh_and_gridtransfer rebuilds the
// mesh but does NOT touch BSSN's deriv workspace (it is solver-owned). Resize it
// to the post-remesh mesh and recompute dt from the new lmax -- mirrors the
// block bssngr_main's GTS driver runs after each remesh. (For the accuracy test
// remesh is OFF, so this never fires; kept for parity + longer smoke runs.)
template <typename TS>
static void bssn_post_remesh_fixup(bssn::BSSNCtx* appCtx, TS* stepper,
                                   int rank_global) {
    bssn::deallocate_bssn_deriv_workspace();
    bssn::allocate_bssn_deriv_workspace(appCtx->get_mesh(), 1);
    stepper->sync_with_mesh();
    appCtx->calculate_full_grid_size();

    ot::Mesh* pmesh = appCtx->get_mesh();
    unsigned int lmin, lmax;
    pmesh->computeMinMaxLevel(lmin, lmax);
    if (!rank_global) printf("New min and max level = (%d, %d)\n", lmin, lmax);

    bssn::BSSN_CURRENT_MIN_DX =
        ((bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
         ((1u << (m_uiMaxDepth - lmax)) / ((double)bssn::BSSN_ELE_ORDER)) /
         ((double)(1u << (m_uiMaxDepth))));
    bssn::BSSN_RK45_TIME_STEP_SIZE = bssn::BSSN_CFL_FACTOR * bssn::BSSN_CURRENT_MIN_DX;

    ts::TSInfo ts_in = appCtx->get_ts_info();
    ts_in._m_uiTh    = bssn::BSSN_RK45_TIME_STEP_SIZE;
    appCtx->set_ts_info(ts_in);
}

int main(int argc, char** argv) {
    // 0- NUTS (LTS)   1- UTS (GTS)
    unsigned int ts_mode = 0;

    if (argc < 2) {
        std::cout << "No parameter file was given, exiting..." << std::endl;
        std::cout << "Usage: " << argv[0] << " paramFile [ts_mode]" << std::endl;
        std::cout << "   ts_mode: 0 = ENUTS/LTS (default), 1 = ETS/GTS control"
                  << std::endl;
        exit(0);
    }

    if (argc > 2) ts_mode = std::atoi(argv[2]);

    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    // 1. read the parameter file.
    if (!rank) std::cout << " reading parameter file :" << argv[1] << std::endl;
    bssn::readParamFile(argv[1], comm);

    // initialize the logger (BSSNCtx::initialize / init_grid log through it).
    dendro::logger::initialize(
        "dendro", bssn::DENDRO_LOG_FILE, 0, bssn::DENDRO_LOG_FILE_LEVEL,
        bssn::DENDRO_LOG_CONSOLE_LEVEL, bssn::DENDRO_LOG_FORCE_FILE_FLUSH);
    dendro::logger::setup_crash_handler();

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

    // 2. generate the initial grid (BSSN puncture initial data). Mirrors
    //    bssngr_main.cpp -- NOT EM4's IC path.
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double, double, double, double*)> f_init =
        [](double x, double y, double z, double* var) {
            bssn::initialDataFunctionWrapper(x, y, z, var);
        };
    std::function<void(double, double, double, double*)> f_init_flat =
        [](double x, double y, double z, double* var) {
            bssn::minkowskiInitialData(x, y, z, var);
        };

    const unsigned int interpVars = bssn::BSSN_NUM_VARS;
    unsigned int varIndex[interpVars];
    for (unsigned int i = 0; i < bssn::BSSN_NUM_VARS; i++) varIndex[i] = i;

    if (bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY) {
        if (!rank)
            std::cout << YLW << "Using block adaptive mesh. AMR disabled " << NRM
                      << std::endl;
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
    bssn::BSSN_CURRENT_MIN_DX =
        ((bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
         ((1u << (m_uiMaxDepth - lmax)) / ((double)bssn::BSSN_ELE_ORDER)) /
         ((double)(1u << (m_uiMaxDepth))));
    bssn::BSSN_RK45_TIME_STEP_SIZE =
        bssn::BSSN_CFL_FACTOR * bssn::BSSN_CURRENT_MIN_DX;
    tmpNodes.clear();

    if (!rank) {
        std::cout << " lmin: " << lmin << " lmax: " << lmax << std::endl;
        std::cout << " ts_mode: " << ts_mode << " (0=ENUTS/LTS, 1=ETS/GTS)"
                  << std::endl;
        std::cout << " dt: " << bssn::BSSN_RK45_TIME_STEP_SIZE << std::endl;
    }

    // RK order: honor BSSN_RK_TYPE, default RK4 (the order this test measures).
    // Both arms use plain RK coefficients (no MSRK) so the GTS control is a
    // clean RK4 self-convergence = 4.00.
    ts::ETSType tsType = ts::ETSType::RK4;
    if ((RKType)bssn::BSSN_RK_TYPE == RKType::RK3)
        tsType = ts::ETSType::RK3;
    else if ((RKType)bssn::BSSN_RK_TYPE == RKType::RK5)
        tsType = ts::ETSType::RK5;

    if (ts_mode == 0) {
        // ================= ENUTS (LTS) =================
        bssn::BSSNCtx* appCtx = new bssn::BSSNCtx(mesh);
        // useWpart=true remeshes at init (weighted LTS repartition) -> ENUTS
        // would solve a DIFFERENT grid than GTS. For a clean LTS-vs-GTS accuracy
        // check, DENDRO_NUTS_NO_WPART=1 keeps ENUTS on the same grid as ETS.
        const bool useWpart = !std::getenv("DENDRO_NUTS_NO_WPART");
        ts::ExplicitNUTS<DendroScalar, bssn::BSSNCtx>* enuts =
            new ts::ExplicitNUTS<DendroScalar, bssn::BSSNCtx>(appCtx, useWpart);

        enuts->set_evolve_vars(appCtx->get_evolution_vars());
        enuts->set_ets_coefficients(tsType);

        const unsigned int rank_global = enuts->get_global_rank();

        for (enuts->init(); enuts->curr_time() < bssn::BSSN_RK_TIME_END;
             enuts->evolve()) {
            const DendroIntL step = enuts->curr_step();

            enuts->dump_load_statistics(std::cout);
            // model the LTS sub-cycling win (throttle-immune: counts, not
            // wall-clock) -- per-level histogram + LTS/GTS work + est. speedup.
            if (step == 0) enuts->dump_est_speedup(std::cout);

            if (!rank_global)
                std::cout << GRN << "[Explicit NUTS]: Executing step :  " << step
                          << std::setw(10)
                          << "\tcurrent time :" << enuts->curr_time()
                          << std::setw(10) << "\t dt(min):" << enuts->get_dt_min()
                          << std::setw(10) << "\t dt(max):" << enuts->get_dt_max()
                          << std::setw(10) << "\t" << NRM << std::endl;

            appCtx->terminal_output();

            bool isRemesh = false;
            if (bssn::BSSN_REMESH_TEST_FREQ > 0 &&
                (step % bssn::BSSN_REMESH_TEST_FREQ) == 0 && step != 0)
                isRemesh = appCtx->is_remesh();

            if (isRemesh) {
                if (!rank_global)
                    std::cout << "[Explicit NUTS]: Remesh triggered" << std::endl;
                appCtx->remesh_and_gridtransfer(bssn::BSSN_DENDRO_GRAIN_SZ,
                                                bssn::BSSN_LOAD_IMB_TOL,
                                                bssn::BSSN_SPLIT_FIX);
                bssn_post_remesh_fixup(appCtx, enuts, rank_global);
            }

            if (bssn::BSSN_IO_OUTPUT_FREQ > 0 &&
                (step % bssn::BSSN_IO_OUTPUT_FREQ) == 0)
                appCtx->write_vtu();

            if (bssn::BSSN_CHECKPT_FREQ > 0 &&
                (step % bssn::BSSN_CHECKPT_FREQ) == 0)
                appCtx->write_checkpt();

            appCtx->prepare_for_next_iter();
        }

        dump_evar_bin(appCtx->get_mesh(), appCtx->get_evolution_vars());

        ot::Mesh* tmp_mesh = appCtx->get_mesh();
        delete appCtx;
        delete tmp_mesh;
        delete enuts;

    } else if (ts_mode == 1) {
        // ================= ETS (GTS, uniform RK) -- the RK4 control =========
        bssn::BSSNCtx* appCtx = new bssn::BSSNCtx(mesh);
        ts::ETS<DendroScalar, bssn::BSSNCtx>* ets =
            new ts::ETS<DendroScalar, bssn::BSSNCtx>(appCtx);
        ets->set_evolve_vars(appCtx->get_evolution_vars());
        ets->set_ets_coefficients(tsType);

        const unsigned int rank_global = ets->get_global_rank();

        for (ets->init(); ets->curr_time() < bssn::BSSN_RK_TIME_END;
             ets->evolve()) {
            const DendroIntL step = ets->curr_step();

            if (!rank_global)
                std::cout << GRN << "[ETS]: Executing step :  " << step
                          << std::setw(10)
                          << "\tcurrent time :" << ets->curr_time()
                          << std::setw(10) << "\t dt:" << ets->ts_size()
                          << std::setw(10) << "\t" << NRM << std::endl;

            appCtx->terminal_output();

            bool isRemesh = false;
            if (bssn::BSSN_REMESH_TEST_FREQ > 0 &&
                (step % bssn::BSSN_REMESH_TEST_FREQ) == 0 && step != 0)
                isRemesh = appCtx->is_remesh();

            if (isRemesh) {
                if (!rank_global)
                    std::cout << "[ETS] : Remesh is triggered.\n";
                appCtx->remesh_and_gridtransfer(bssn::BSSN_DENDRO_GRAIN_SZ,
                                                bssn::BSSN_LOAD_IMB_TOL,
                                                bssn::BSSN_SPLIT_FIX);
                bssn_post_remesh_fixup(appCtx, ets, rank_global);
            }

            if (bssn::BSSN_IO_OUTPUT_FREQ > 0 &&
                (step % bssn::BSSN_IO_OUTPUT_FREQ) == 0)
                appCtx->write_vtu();

            if (bssn::BSSN_CHECKPT_FREQ > 0 &&
                (step % bssn::BSSN_CHECKPT_FREQ) == 0)
                appCtx->write_checkpt();

            appCtx->prepare_for_next_iter();
        }

        dump_evar_bin(appCtx->get_mesh(), appCtx->get_evolution_vars());

        ot::Mesh* tmp_mesh = appCtx->get_mesh();
        delete appCtx;
        delete tmp_mesh;
        delete ets;
    }

    MPI_Finalize();

    if (!rank) std::cout << GRN << "bssnSolverNUTS finished!" << NRM << std::endl;

    return 0;
}

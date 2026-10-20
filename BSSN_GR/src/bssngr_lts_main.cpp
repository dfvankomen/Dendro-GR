/**
 * @file bssngr_lts_main.cpp
 * @brief bssngr_main.cpp with the local time stepper as a run-time choice.
 *
 * TS_MODE 1 (default) is the uniform stepper, step for step what
 * bssngr_main.cpp does. TS_MODE 0 builds ts::ExplicitNUTS on the same context
 * and runs the same loop: tracker, merger switch, remesh, extraction, output.
 * ExplicitNUTS derives from ts::ETS, so only the construction differs; the
 * step-count frequencies become simulation-time intervals under LTS because
 * one LTS step is a coarse step.
 */
#include <algorithm>
#include <chrono>
#include <cstdlib>  // std::getenv for the OMP_SCHEDULE default guard
#include <ctime>  // time_t / time() / ctime() — gcc 15 + <chrono> doesn't drag these in
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <vector>
#if defined(DENDRO_HYBRID_OMP) && defined(_OPENMP)
#include <omp.h>
#endif

#include "TreeNode.h"
#include "aeh.h"
#include "bssnAEH.h"
#include "bssnCtx.h"
#include "gr.h"
#include "grUtils.h"
#include "logger.h"
#include "mathUtils.h"  // normL2
#include "mesh.h"
#include "meshUtils.h"
#include "mpi.h"
#include "octUtils.h"
#include "parameters.h"
#include "sdc.h"

// DENDRO_EVAR_DUMP=<prefix>: write this rank's local evolution vector to
// <prefix>_rank<globalRank>.bin at the end of the run (pointwise comparison of
// two runs on the same fixed mesh).
static void dump_evar_bin(const ot::Mesh* mesh, bssn::DVec& evar) {
    const char* pfx = std::getenv("DENDRO_EVAR_DUMP");
    if (!pfx || !mesh->isActive()) return;
    if (std::getenv("DENDRO_EVAR_DUMP_POS")) {
        // <prefix>_pos_rank<r>.bin: one record per local CG node, in the order
        // of the plain dump: xyz (3 x uint64, octree coordinate times the
        // element order), dup flag and owner-element level (2 x uint32).
        // Localises a dump difference by position and level.
        const auto& c2d = mesh->getCG2DGMap();
        const auto& el  = mesh->getAllElements();
        const auto& dup = mesh->getLocalNonPrimaryNodes();
        const unsigned int eo = mesh->getElementOrder(), md = m_uiMaxDepth;
        char fnp[512];
        snprintf(fnp, sizeof(fnp), "%s_pos_rank%d.bin", pfx,
                 mesh->getMPIRankGlobal());
        FILE* fp = fopen(fnp, "wb");
        for (unsigned int cg = mesh->getNodeLocalBegin();
             fp && cg < mesh->getNodeLocalEnd(); cg++) {
            unsigned int e, i, j, k;
            mesh->dg2eijk(c2d[cg], e, i, j, k);
            const unsigned long long len = 1ull << (md - el[e].getLevel());
            unsigned long long xyz[3] = {
                el[e].getX() * (unsigned long long)eo + i * len,
                el[e].getY() * (unsigned long long)eo + j * len,
                el[e].getZ() * (unsigned long long)eo + k * len};
            unsigned int meta[2] = {
                (unsigned)std::binary_search(dup.begin(), dup.end(), cg),
                el[e].getLevel()};
            fwrite(xyz, 8, 3, fp);
            fwrite(meta, 4, 2, fp);
        }
        if (fp) fclose(fp);
    }
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

// BSSN_CONSTRAINT_MONITOR_DT = t: print the constraint L2 norms at the first
// step at or past every multiple of t (0 = off).
static void dump_constraint_norms(bssn::BSSNCtx* ctx, DendroIntL step, double t,
                                  const char* tag) {
    static double next = 0;
    const double dt    = bssn::BSSN_CONSTRAINT_MONITOR_DT;
    if (dt <= 0 || t < next - 1e-9 * std::max(1.0, t)) return;
    next        = (std::floor(t / dt + 1e-9) + 1.0) * dt;

    ot::Mesh* m = ctx->get_mesh();
    if (!m->isActive()) return;
    ctx->compute_constraint_variables();
    bssn::DVec& cv = ctx->get_constraint_vars();
    DendroScalar* c[bssn::BSSN_CONSTRAINT_NUM_VARS];
    cv.to_2d(c);
    const unsigned int nlb  = m->getNodeLocalBegin();
    const unsigned int nloc = m->getNumLocalMeshNodes();
    MPI_Comm comm           = m->getMPICommunicator();
    DendroScalar h  = normL2(c[bssn::VAR_CONSTRAINT::C_HAM] + nlb, nloc, comm);
    DendroScalar m0 = normL2(c[bssn::VAR_CONSTRAINT::C_MOM0] + nlb, nloc, comm);
    DendroScalar m1 = normL2(c[bssn::VAR_CONSTRAINT::C_MOM1] + nlb, nloc, comm);
    DendroScalar m2 = normL2(c[bssn::VAR_CONSTRAINT::C_MOM2] + nlb, nloc, comm);
    DendroScalar mom = std::sqrt(m0 * m0 + m1 * m1 + m2 * m2);
    if (!m->getMPIRankGlobal())
        printf("[CONSTRAINT] %s step %lld t %.6f ||H||=%.8e ||M||=%.8e\n", tag,
               (long long)step, t, (double)h, (double)mom);
}

int main(int argc, char** argv) {
    // 0- NUTS 1-UTS
    unsigned int ts_mode = 1;

    if (argc < 2) {
        std::cout << "usage: " << argv[0] << " PARAM_FILE [TS_MODE]"
                  << std::endl;
        std::cout << std::endl << "options:" << std::endl;
        std::cout << "  PARAM_FILE" << std::endl
                  << "      Path to the parameter file (.toml)" << std::endl;
        std::cout << "  TS_MODE" << std::endl
                  << "      Time stepper mode." << std::endl;
        std::cout << "        0 - local time stepping (LTS)" << std::endl;
        std::cout << "        1 - global time stepping (" << GRN << "default"
                  << NRM << ")" << std::endl;
        return 0;
    }

    if (argc > 2) ts_mode = std::atoi(argv[2]);

    // Hybrid OpenMP/MPI needs at least FUNNELED: only the master thread calls
    // MPI (threading happens inside compute regions, not around MPI calls).
    // Requested unconditionally -- harmless for pure-MPI runs.
    int mpi_required = MPI_THREAD_FUNNELED;
    int mpi_provided = 0;
    MPI_Init_thread(&argc, &argv, mpi_required, &mpi_provided);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    if (!rank && mpi_provided < mpi_required) {
        std::cout << YLW
                  << "WARNING: MPI provided thread level " << mpi_provided
                  << " < requested MPI_THREAD_FUNNELED; the hybrid OpenMP/MPI "
                     "path may be unsafe with this MPI build."
                  << NRM << std::endl;
    }

#if defined(DENDRO_HYBRID_OMP) && defined(_OPENMP)
    // Early default for the RHS loop's schedule(runtime): dynamic,1 (historical)
    // when OMP_SCHEDULE is unset. The authoritative resolution (which also honors
    // the parfile BSSN_HYBRID_RHS_SCHEDULE) happens after readParamFile via
    // bssn::set_rhs_omp_schedule(); the effective value is logged there.
    if (std::getenv("OMP_SCHEDULE") == nullptr)
        omp_set_schedule(omp_sched_dynamic, 1);
#endif

    // --- CPU capability gate -----------------------------------------
    // If this binary was compiled for AVX2/AVX-512 but the runtime CPU
    // lacks those features (e.g., built on a login node, run on an
    // older compute node), fail fast with a clear message instead of
    // crashing on illegal-instruction. __builtin_cpu_supports requires
    // a string literal so we check each feature by name directly.
#if defined(__GNUC__) || defined(__clang__)
    __builtin_cpu_init();
    {
        const char* missing = nullptr;
 #if defined(BSSN_USE_CASCADE_AVX) || defined(BSSN_USE_CASCADE_AVX_FUSED) \
  || defined(BSSN_USE_CASCADE_AVX512) || defined(BSSN_USE_CASCADE_AVX512_FUSED)
        if (!__builtin_cpu_supports("avx2")) missing = "avx2";
        else if (!__builtin_cpu_supports("fma")) missing = "fma";
 #endif
 #if defined(BSSN_USE_CASCADE_AVX512) || defined(BSSN_USE_CASCADE_AVX512_FUSED)
        if (!missing && !__builtin_cpu_supports("avx512f"))  missing = "avx512f";
        if (!missing && !__builtin_cpu_supports("avx512dq")) missing = "avx512dq";
 #endif
        if (missing) {
            if (!rank) {
                std::cerr << "FATAL: binary requires CPU feature '" << missing
                          << "' but runtime CPU lacks it.\n"
                          << "       Rebuild for the target node's CPU, or "
                             "run on a node that supports " << missing << ".\n";
            }
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
#endif
    // ------------------------------------------------------------------

    if (!rank) {
        std::cout << "======================================" << std::endl;
        std::cout << GRN << ":::: Now initializing BSSN Solver ::::" << NRM
                  << std::endl;
        if (ts_mode == 0) {
            std::cout << YLW
                      << "      - Running with the Non-Uniform/Spatially "
                         "Adaptive Time Stepper (NUTS/SATS)"
                      << NRM << std::endl;
        } else {
            std::cout << YLW
                      << "      - Running with the Uniform Time Stepper (UTS)"
                      << NRM << std::endl;
        }

        std::cout << "======================================" << std::endl;
    }

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

    std::vector<std::string> arg_s(argv, argv + argc);
    bssn::printGitInformation(rank, arg_s);

    // 1 . read the parameter file.
    if (!rank) std::cout << " reading parameter file :" << argv[1] << std::endl;
    bssn::readParamFile(argv[1], comm);
    bssn::resolve_cahd_form(ts_mode == 0, comm);

    if (ts_mode == 0) {
        // LTS defaults that the stepper, the step rule and the block loop read
        // on first use, so they are set before anything is built. each one
        // stays as the environment gives it if set.
        // owner sync on the level-wise partition
        setenv("DENDRO_ENUTS_OWNER_SYNC_LABELLED", "1", 0);
        // a rank with fewer than 16 blocks on its finest level splits them, so
        // the finest-only sub-steps keep its threads busy (level-wise mesh only;
        // 0 turns it off)
        setenv("DENDRO_OCT2BLK_MIN_BLOCKS", "16", 0);
#ifdef BSSN_WITH_GRAPH_PARTITIONING
        // a par file that does not name a partition method gets the level-wise
        // one under LTS; an explicit value (0 = SFC included) is kept
        {
            int named = 0;
            if (!rank) {
                std::ifstream pf(argv[1]);
                std::string line;
                const std::regex key("^\\s*BSSN_PARTITIONING_METHOD\\s*=");
                while (std::getline(pf, line))
                    if (std::regex_search(line, key)) named = 1;
            }
            MPI_Bcast(&named, 1, MPI_INT, 0, comm);
            if (!named && bssn::BSSN_PARTITIONING_METHOD == 0) {
                bssn::BSSN_PARTITIONING_METHOD = 4;
                if (!rank)
                    std::cout << YLW
                              << "[LTS] par file names no partition method: "
                                 "using the level-wise partition "
                                 "(BSSN_PARTITIONING_METHOD = 4)"
                              << NRM << std::endl;
            }
        }
#endif
    }

#if defined(DENDRO_HYBRID_OMP) && defined(_OPENMP)
    // Resolve the RHS block-loop schedule now that the parfile is parsed:
    // BSSN_HYBRID_RHS_SCHEDULE (parfile) > OMP_SCHEDULE (env) > dynamic,1.
    // (Also re-applied in allocate_bssn_deriv_workspace(); this is the early,
    // authoritative resolution + the startup log line for A/B runs.)
    bssn::set_rhs_omp_schedule();
    if (!rank) {
        omp_sched_t sched_kind;
        int sched_chunk;
        omp_get_schedule(&sched_kind, &sched_chunk);
        // Mask the omp_sched_monotonic modifier (high bit) that OMP_SCHEDULE=...
        // sets, so the base kind reads cleanly (1=static 2=dynamic 3=guided
        // 4=auto).
        const int base = (int)sched_kind & 0x7fffffff;
        const char* name = (base == 1)   ? "static"
                           : (base == 2) ? "dynamic"
                           : (base == 3) ? "guided"
                           : (base == 4) ? "auto"
                                         : "unknown";
        std::cout << "[hybrid] RHS block schedule (effective): " << name
                  << " chunk=" << sched_chunk << "  (kind_raw=" << (int)sched_kind
                  << "; set via parfile BSSN_HYBRID_RHS_SCHEDULE or OMP_SCHEDULE; "
                     "bare 'static' = NUMA-local contiguous partition -- see the "
                     "RHS NUMA Tax)"
                  << std::endl;
    }
#endif

    // get the current time
    auto now    = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);

    std::ostringstream tmp_stringstream;
    tmp_stringstream << std::put_time(std::localtime(&time_t),
                                      "%Y-%m-%d-%H-%M-%S");
    std::string timestamp = tmp_stringstream.str();

    // then write it right back out:
    std::string outputFileName =
        bssn::BSSN_PROFILE_FILE_PREFIX + "__PARAM_DUMP__" + timestamp + ".toml";
    bssn::writeParamTOMLFile(outputFileName.c_str(), comm);

    // initialize the logger
    dendro::logger::initialize(
        "dendro", bssn::DENDRO_LOG_FILE, 0, bssn::DENDRO_LOG_FILE_LEVEL,
        bssn::DENDRO_LOG_CONSOLE_LEVEL, bssn::DENDRO_LOG_FORCE_FILE_FLUSH);
    // this call sets up the ability for the dendro logger to force a flush when
    // a signal is received. If, for some reason, this appears to be "not
    // enough" for serious debugging, make sure to set FORCE_FILE_FLUSH which
    // will force writing to file each time we make a debug log
    dendro::logger::setup_crash_handler();
    dendro::logger::info(
        "Logger initialized, now initializing the rest of the data structures");
    dendro::logger::info("All parameters are also read!");

    // barrier on the dumping of parameter file to avoid gumming up the console
    // output, it'll take just a second
    MPI_Barrier(comm);
    int root = std::min(1, npes - 1);
    bssn::dumpParamFile(std::cout, root, comm);
    MPI_Barrier(comm);

    _InitializeHcurve(bssn::BSSN_DIM);
    m_uiMaxDepth = bssn::BSSN_MAXDEPTH;

    if (bssn::BSSN_NUM_VARS % bssn::BSSN_ASYNC_COMM_K != 0) {
        if (!rank)
            std::cout << "[overlap communication error]: total BSSN_NUM_VARS: "
                      << bssn::BSSN_NUM_VARS
                      << " is not divisable by BSSN_ASYNC_COMM_K: "
                      << bssn::BSSN_ASYNC_COMM_K << std::endl;
        MPI_Abort(comm, 1);
    }

    // NOTE: this is where we originally had the check that Extract freq needed
    // to be greater than IO freq, but decoupling the parameters means that this
    // isn't "required", even though it's probably recommended.

    // 2. generate the initial grid.
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double, double, double, double*)> f_init =
        [](double x, double y, double z, double* var) {
            bssn::initialDataFunctionWrapper(x, y, z, var);
        };
    std::function<double(double, double, double)> f_init_alpha =
        [](double x, double y, double z) {
            double var[24];
            bssn::initialDataFunctionWrapper(x, y, z, var);
            return var[0];
        };

    const unsigned int interpVars = bssn::BSSN_NUM_VARS;
    unsigned int varIndex[interpVars];
    for (unsigned int i = 0; i < bssn::BSSN_NUM_VARS; i++) varIndex[i] = i;

    /*varIndex[0]=bssn::VAR::U_ALPHA;
    varIndex[1]=bssn::VAR::U_CHI;*/
    DendroIntL localSz, globalSz;
    double t_stat;
    double t_stat_g[3];

    bssn::timer::t_f2o.start();

    if (bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY) {
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

            MPI_Abort(comm, 1);
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

    if (ts_mode == 0 || ts_mode == 1) {
        const bool useLTS = (ts_mode == 0);
        if (!rank)
            std::cout << GRN << "Now setting up the "
                      << (useLTS ? "local" : "uniform") << " time stepper!"
                      << NRM << std::endl;
        bssn::BSSNCtx::set_subcycled(useLTS);
        bssn::BSSNCtx* bssnCtx = new bssn::BSSNCtx(mesh);

        // select time integrator based on BSSN_RK_TYPE (index-matched to
        // ts::ETSType):
        //   0=RK3, 1=RK4, 2=RK5 (Butcher 6-stage)
        //   3=MSRK2_1, 4=MSRK2_2, 5=MSRK3 (multistep RK, arXiv:2603.05763)
        //   6=RALSTON, 7=CASH_KARP, 8=RKF45, 9=NYSTROM, 10=RK6 (single-step)
        //   11=RK6_TSRK (general multistep RK, order 6, depth 2)
        // MSRK/TSRK methods reuse RHS evals from previous steps to cut ghost
        // exchanges (~20-40% speedup); TSRK stores 8 extra derivative clouds.
        const RKType rkType = (RKType)bssn::BSSN_RK_TYPE;

        ts::ETS<DendroScalar, bssn::BSSNCtx>* ets = nullptr;
        // non-null only under LTS, for the calls ETS does not have
        ts::ExplicitNUTS<DendroScalar, bssn::BSSNCtx>* enuts = nullptr;

        if (useLTS && bssn::BSSN_PARTITIONING_METHOD != 0) {
            // the labelling replaces the weighted start-up partition
            setenv("DENDRO_NUTS_NO_WPART", "1", 1);
        }
        if (useLTS) {
            // LTS equals the uniform stepper at step ratio 1 only with the
            // interface fixes and the per-stage time, so they are on here
            // unless the environment says otherwise.
            setenv("DENDRO_ENUTS_OWNER_SYNC", "1", 0);
            setenv("DENDRO_ENUTS_SYNC_FAST", "1", 0);
            setenv("DENDRO_ENUTS_STAGE_TIME", "1", 0);
            // owner-sync cuts measured on q12 (loop 61.6 -> 54.7 s): no interior
            // self-copy, no per-sub-step element-time exchange, threaded pack /
            // unpack, and one shared parent-to-child store per pass
            setenv("DENDRO_ENUTS_OWN_NO_SELFCOPY", "1", 0);
            setenv("DENDRO_ENUTS_OWN_NO_ETIME_XCH", "1", 0);
            setenv("DENDRO_ENUTS_OWN_XCH_OMP", "1", 0);
            setenv("DENDRO_ENUTS_OWN_P2C_SHARED", "2", 0);
            // send a ghost slot only when a block that reads it steps (q12: 58 %
            // fewer doubles exchanged, loop 50.9 -> 46.8 s)
            setenv("DENDRO_ENUTS_XCH_SCOPED", "1", 0);
            // point location for interpolation (tracker, wave extraction) by
            // containment only: the splitter route assumes one SFC range per
            // rank and misses on a labelled mesh (same values, no misses)
            setenv("DENDRO_INTERP_DIRECT", "1", 0);
            // a step is a coarse step here and the finder is gated on time,
            // so its own step-count gate on the diagnostics files is wrong
            if (AEH::ah_bah) AEH::ah_bah->set_write_every_find(true);
            // the work-weighted start-up partition: apply it after the initial
            // grid has converged (so weighted and unweighted runs start from one
            // mesh) and with the evolution vector's ghosts filled first
            setenv("DENDRO_ENUTS_WPART_LATE", "1", 0);
            setenv("DENDRO_ENUTS_WPART_GHOST_SYNC", "1", 0);
            // the weighted repartition remeshes at init; DENDRO_NUTS_NO_WPART=1
            // keeps the LTS run on the grid the uniform run starts from
            const bool useWpart = !std::getenv("DENDRO_NUTS_NO_WPART");
            ts::LtsOptions opts;
            opts.iface      = bssn::BSSN_LTS_IFACE;
            opts.g_exchange = bssn::BSSN_LTS_G_EXCHANGE;
            opts.ratio_cap  = bssn::BSSN_LTS_RATIO_CAP;
            enuts = new ts::ExplicitNUTS<DendroScalar, bssn::BSSNCtx>(
                bssnCtx, useWpart, opts);
            ets   = enuts;
            ets->set_evolve_vars(bssnCtx->get_evolution_vars());
            // RK5 segfaults in the first LTS step (any partition, any switches;
            // gate seq36, 2026-10-14), so it is refused until the stepper is fixed
            if (rkType != RKType::RK3 && rkType != RKType::RK4) {
                if (!rank)
                    std::cerr << "LTS supports RK3 / RK4, not BSSN_RK_TYPE="
                              << bssn::BSSN_RK_TYPE << std::endl;
                MPI_Abort(comm, 1);
            }
            ets->set_ets_coefficients((ts::ETSType)rkType);
        } else if (rkType == RKType::RK4_MSRK2_1 || rkType == RKType::RK4_MSRK2_2 ||
            rkType == RKType::RK4_MSRK3) {
            ts::ETSType msrkVariant;
            if (rkType == RKType::RK4_MSRK2_1)
                msrkVariant = ts::ETSType::RK4_MSRK2_1;
            else if (rkType == RKType::RK4_MSRK2_2)
                msrkVariant = ts::ETSType::RK4_MSRK2_2;
            else
                msrkVariant = ts::ETSType::RK4_MSRK3;

            ets = new ts::ETS_MSRK<DendroScalar, bssn::BSSNCtx>(
                bssnCtx, msrkVariant);
        } else if (rkType == RKType::RK6_TSRK) {
            // TSRK takes no variant arg: fixed to the order-6 tableau, and
            // self-bootstraps its history via the base RK6 stepper.
            ets = new ts::ETS_TSRK<DendroScalar, bssn::BSSNCtx>(bssnCtx);
        } else {
            ets = new ts::ETS<DendroScalar, bssn::BSSNCtx>(bssnCtx);

            // RKType is index-matched to ts::ETSType, so a direct cast covers
            // every single-step method (RK3/4/5 + RALSTON/CASH_KARP/RKF45/
            // NYSTROM/RK6).
            if (ets->set_ets_coefficients((ts::ETSType)rkType) != 0) {
                if (!rank)
                    std::cerr << "Unknown BSSN_RK_TYPE=" << bssn::BSSN_RK_TYPE
                              << ", defaulting to RK4" << std::endl;
                ets->set_ets_coefficients(ts::ETSType::RK4);
            }
        }

        ets->set_evolve_vars(bssnCtx->get_evolution_vars());
        ets->init();
#ifdef BSSN_WITH_GRAPH_PARTITIONING
        // repartition once the initial data is set: build a relabelled twin of
        // the SFC mesh, move the evolution vector onto it, re-sync the stepper
        if (bssn::BSSN_PARTITIONING_METHOD != 0) {
            if (!rank)
                std::cout << CYN << BLD << "Applying partitioning method "
                          << bssn::BSSN_PARTITIONING_METHOD << NRM << std::endl;
            if (ts::apply_partitioning(*bssnCtx, *ets) && !rank)
                std::cout << GRN << BLD << "...labelled partition active."
                          << NRM << std::endl;
        }
#else
        if (bssn::BSSN_PARTITIONING_METHOD != 0) {
            if (!rank)
                std::cerr << "BSSN_PARTITIONING_METHOD needs a build with "
                             "BSSN_WITH_GRAPH_PARTITIONING" << std::endl;
            MPI_Abort(comm, 1);
        }
#endif
#if defined __PROFILE_CTX__ && defined __PROFILE_ETS__
        std::ofstream outfile;
        char fname[256];
        sprintf(fname, "bssnCtx_%d.txt", npes);
        if (!rank) {
            outfile.open(fname, std::ios_base::app);
            std::time_t now = std::time(nullptr);
            // convert now to string form
            char* dt        = std::ctime(&now);
            outfile << "======================================================="
                       "====="
                    << std::endl;
            outfile << "Current time : " << dt << " --- " << std::endl;
            outfile << "======================================================="
                       "====="
                    << std::endl;
        }

        ets->init_pt();
        bssnCtx->reset_pt();
        ets->dump_pt(outfile);
        // bssnCtx->dump_pt(outfile);
#endif
        double t1              = MPI_Wtime();

        bssnaeh::initialize_aeh();

        // capture the curr step
        const DendroIntL start_step = ets->curr_step();

        ts::Cadence cRemesh, cOut, cGW, cIO, cAEH, cChk;
        // the remesh test skips step 0, so its first interval starts there
        cRemesh.next = bssn::BSSN_REMESH_TEST_FREQ * bssn::BSSN_RK45_TIME_STEP_SIZE;
        auto due = [&](ts::Cadence& c, unsigned int freq, DendroIntL step,
                       double time) {
            return c.due(freq, step, time, bssn::BSSN_RK45_TIME_STEP_SIZE,
                         useLTS);
        };
        if (useLTS && start_step != 0) {
            // restart: an uninterrupted run fires at this step only for an
            // interval boundary inside the coarse step that just ended, so
            // start every gate from the first boundary after its beginning
            const double tPrev = ets->curr_time() - enuts->get_dt_max();
            // a run restored after the merger was already on the post-merger
            // intervals (the loop sets them again at its first step)
            if (bssnCtx->is_bh_merged(0.1)) {
                bssn::BSSN_REMESH_TEST_FREQ =
                    bssn::BSSN_REMESH_TEST_FREQ_AFTER_MERGER;
                bssn::BSSN_GW_EXTRACT_FREQ =
                    bssn::BSSN_GW_EXTRACT_FREQ_AFTER_MERGER;
            }
            const double dtf = bssn::BSSN_RK45_TIME_STEP_SIZE;
            cRemesh.resume(bssn::BSSN_REMESH_TEST_FREQ, tPrev, dtf);
            cOut.resume(bssn::BSSN_TIME_STEP_OUTPUT_FREQ, tPrev, dtf);
            cGW.resume(bssn::BSSN_GW_EXTRACT_FREQ_TRUE, tPrev, dtf);
            cIO.resume(bssn::BSSN_IO_OUTPUT_FREQ_TRUE, tPrev, dtf);
            cAEH.resume(AEH::AEH_SOLVER_FREQ, tPrev, dtf);
            cChk.resume(bssn::BSSN_CHECKPT_FREQ, tPrev, dtf);
        }

        dendro::logger::info("Now beginning ETS time stepper");
        while (ets->curr_time() < bssn::BSSN_RK_TIME_END) {
            const DendroIntL step            = ets->curr_step();
            const DendroScalar time          = ets->curr_time();

            bssn::BSSN_CURRENT_RK_COORD_TIME = time;
            bssn::BSSN_CURRENT_RK_STEP       = step;

            const bool isActive              = ets->is_active();
            const unsigned int rank_global   = ets->get_global_rank();

            // things that should happen **only** on time step 0
            if (step == 0) {
                if (!rank_global) {
                    std::cout << BLU
                              << "[ETS] : Timestep 0 - ensuring a few things "
                                 "are taken care of..."
                              << NRM << std::endl;
                }
                // for our scaling operation, we want to make sure that the
                // constraints are computed and handled
                bssnCtx->compute_constraint_variables();

                // grid summary .dat always written (not profiling telemetry)
                bssnCtx->write_grid_summary_data();

                if (!rank_global) {
                    std::cout << BLU
                              << "[ETS] : Timestep 0 - Finished with things "
                                 "that should always be done at time 0!"
                              << NRM << std::endl;
                }
            }

            // on restore, but only on restore and not at the beginning
            if (step != 0 && step == start_step) {
                if (!rank_global)
                    std::cout << BLD << GRN
                              << "[ETS] : CHECKPOINT RESTORED. Doing "
                                 "additional cleanup.\n"
                              << NRM << std::endl;
                // there's no guarantee that the constraints will be computed
                // right on restore, due to the timing. we still need them
                // populated with "good" data
                bssnCtx->compute_constraint_variables();
            }

            const bool is_merged = bssnCtx->is_bh_merged(0.1);
            if (is_merged) {
                // make sure we set that bh is merged!
                // NOTE: don't worry, this won't update the bssn ctx object
                // after the first time, unless something modifies another
                // internal variable
                bssnCtx->set_is_merged(time, step);

                // bssn::BSSN_REMESH_TEST_FREQ=3 *
                // bssn::BSSN_REMESH_TEST_FREQ_AFTER_MERGER;
                // bssn::BSSN_MINDEPTH=5;
                // TODO: make BSSN refinement mode POST MERGER an option!

                // wkb 5 Sept 2024: disable these two lines
                // so I can test other refinement modes
                // bssn::BSSN_REFINEMENT_MODE = bssn::RefinementMode::WAMR;
                // bssn::BSSN_USE_WAVELET_TOL_FUNCTION = 1;
                bssn::BSSN_REMESH_TEST_FREQ =
                    bssn::BSSN_REMESH_TEST_FREQ_AFTER_MERGER;
                bssn::BSSN_GW_EXTRACT_FREQ =
                    bssn::BSSN_GW_EXTRACT_FREQ_AFTER_MERGER;

                // ONLY ENABLE CAKO DURING MERGER
                if (bssn::BSSN_KO_SIGMA_SCALE_BY_CONFORMAL_POST_MERGER_ONLY) {
                    bssn::BSSN_CAKO_ENABLED = true;
                }
            }

            // the uniform stepper remeshes with puncture positions one fine
            // step old; under LTS they would be a whole coarse step old, so
            // bring them up to date first (the later call is then a no-op)
            if (useLTS) bssnCtx->evolve_bh_loc();

            if (step != 0 && due(cRemesh, bssn::BSSN_REMESH_TEST_FREQ, step, time)) {
                ts::remesh_if_needed(*bssnCtx, *ets, bssn::BSSN_DENDRO_GRAIN_SZ,
                                     bssn::BSSN_LOAD_IMB_TOL, bssn::BSSN_SPLIT_FIX);
                // grid summary .dat always written, even after a remesh
                bssnCtx->write_grid_summary_data();
            }

            // print terminal output
            if (due(cOut, bssn::BSSN_TIME_STEP_OUTPUT_FREQ, step, time)) {
                if (!rank_global)
                    std::cout << BLD << GRN << "[ETS - BSSN] : SOLVER UPDATE\n"
                              << NRM << "\tCurrent Step: " << ets->curr_step()
                              << "\t\tCurrent time: " << ets->curr_time()
                              << "\tdt: " << ets->ts_size() << "\t"
                              << std::endl;

                if (enuts && !rank_global)
                    std::cout << "\t[LTS] dt(min): " << enuts->get_dt_min()
                              << "\tdt(max): " << enuts->get_dt_max()
                              << std::endl;

                bssnCtx->terminal_output();
            }
            if (enuts && step == 0) enuts->dump_est_speedup(std::cout);
            dump_constraint_norms(bssnCtx, step, time, useLTS ? "LTS" : "GTS");

            // Per-evolve profile snapshot, once per evolve() interval (step %
            // FREQ misfires since evolve() advances several sub-steps per call).
#ifdef ENABLE_DENDRO_PROFILE_COUNTERS
            // Skip the first iteration (no RHS run yet -> all-zero record) but
            // still flush below. emit_profile is rank-identical, so the
            // collective reductions in profileInfo* stay all-or-none.
            const bool emit_profile = (step > start_step);
            if (emit_profile)
                bssn::timer::profileInfoIntermediate(
                    bssn::BSSN_PROFILE_FILE_PREFIX.c_str(), ets->get_mesh(),
                    step);
#if defined(__PROFILE_ETS__) && defined(__PROFILE_CTX__)
            if (emit_profile)
                bssn::timer::profileInfoJSON(
                    bssn::BSSN_PROFILE_FILE_PREFIX.c_str(), ets->get_mesh(), step,
                    &ets->m_uiCtxpt, &bssnCtx->m_uiCtxpt);
            ets->reset_pt();
#else
            if (emit_profile)
                bssn::timer::profileInfoJSON(
                    bssn::BSSN_PROFILE_FILE_PREFIX.c_str(), ets->get_mesh(),
                    step);
#endif
            bssn::timer::resetSnapshot();
#endif  // ENABLE_DENDRO_PROFILE_COUNTERS

            // wkb: update BH locations always
            bssnCtx->evolve_bh_loc();

            if (due(cGW, bssn::BSSN_GW_EXTRACT_FREQ_TRUE, step, time)) {
                if (!rank_global)
                    std::cout << "    Now extracting constraints and GW."
                              << std::endl;

                // evolving the black holes always stores the updated
                // information
                bssnCtx->extract_constraints();
                bssnCtx->extract_gravitational_waves();
                // write bh coordinates at GW extraction time
                bssnCtx->write_bh_coords();
            }

            // Write VTU and BHLocation files
            if (due(cIO, bssn::BSSN_IO_OUTPUT_FREQ_TRUE, step, time)) {
                // this is all IO output, except for extracting the GW waves,
                // which are "independent"

                // write to vtu, which includes writing the BH location data
                bssnCtx->write_vtu();
            }

            // Run AH solver
            if (due(cAEH, AEH::AEH_SOLVER_FREQ, step, time)) {
                // bssnaeh::perform_aeh_step(bssnCtx, rank);

                bssnCtx->findAH();
            }

            ets->evolve();

            // Write checkpoint  data
            if (due(cChk, bssn::BSSN_CHECKPT_FREQ, step, time)) {
                // the block-cadence puncture tracks are not in the checkpoint:
                // fold them in now, so the stored positions are at this time
                // (the call at the top of the next step, or the first one
                // after a restore, then finds nothing left to cover)
                if (useLTS) {
                    bssnCtx->prepare_for_next_iter();
                    bssnCtx->evolve_bh_loc();
                }
                // a step is a coarse step under LTS, so the step count does
                // not alternate the two slots: use the interval index
                int slot = -1;
                if (useLTS)
                    slot = (int)((long long)std::floor(
                                     time / (bssn::BSSN_CHECKPT_FREQ *
                                             bssn::BSSN_RK45_TIME_STEP_SIZE) +
                                     1e-9) %
                                 2);
                ts::checkpoint(*bssnCtx, *ets,
                               [&] { bssnCtx->write_checkpt(slot); });
            }

            bssnCtx->prepare_for_next_iter();
        }

#if defined __PROFILE_CTX__ && defined __PROFILE_ETS__
        ets->dump_pt(outfile);
        // bssnCtx->dump_pt(outfile);
#endif

        dump_evar_bin(bssnCtx->get_mesh(), bssnCtx->get_evolution_vars());
        bssnCtx->lts_rank_work_report();

        double t2 = MPI_Wtime() - t1;
        double t2_g;
        par::Mpi_Allreduce(&t2, &t2_g, 1, MPI_MAX, ets->get_global_comm());
        if (!(ets->get_global_rank()))
            std::cout << " ETS time (max) : " << t2_g << std::endl;

        // Phase-level timing breakdown (max over ranks, accumulated over run)
        {
            auto max_over_ranks = [&](double v) {
                double g;
                par::Mpi_Allreduce(&v, &g, 1, MPI_MAX, ets->get_global_comm());
                return g;
            };
            const double t_rhs_g          = max_over_ranks(bssn::timer::t_rhs.seconds);
            const double t_deriv_g        = max_over_ranks(bssn::timer::t_deriv.seconds);
            const double t_bdyc_g         = max_over_ranks(bssn::timer::t_bdyc.seconds);
            const double t_rkStep_g       = max_over_ranks(bssn::timer::t_rkStep.seconds);
            const double t_unzip_sync_g   = max_over_ranks(bssn::timer::t_unzip_sync.seconds);
            const double t_zip_g          = max_over_ranks(bssn::timer::t_zip.seconds);
            const double t_ioVtu_g        = max_over_ranks(bssn::timer::t_ioVtu.seconds);
            if (!(ets->get_global_rank())) {
                // unzip_sync includes ghost-exchange (bundled in Ctx::unzip).
                // Per-step + UNZIP_WCOMM/UNZIP split available in the JSONL.
                std::cout << " ---- phase breakdown (max over ranks, s) ----" << std::endl;
                std::cout << "   rkStep       : " << t_rkStep_g << std::endl;
                std::cout << "     deriv      : " << t_deriv_g << std::endl;
                std::cout << "     rhs        : " << t_rhs_g << std::endl;
                std::cout << "     bdyc       : " << t_bdyc_g << std::endl;
                std::cout << "   unzip_sync   : " << t_unzip_sync_g << std::endl;
                std::cout << "   zip          : " << t_zip_g << std::endl;
                std::cout << "   ioVtu        : " << t_ioVtu_g << std::endl;
            }
        }

        // BSSNCtx::~BSSNCtx() dereferences m_uiMesh, so mesh must outlive it.
        ot::Mesh* mesh = bssnCtx->get_mesh();
        delete ets;
        delete bssnCtx;
        delete mesh;

    } else {
        std::cout << RED << "Not starting solver, ts_mode needs to be 0 (LTS) or 1!"
                  << NRM << std::endl;
    }

    MPI_Finalize();
    return 0;
}

/**
 * @file bssn_scaling_bench.cpp
 * @brief Strong / weak scaling benchmark for the BSSN hybrid (MPI+OpenMP) path.
 *
 * Builds either a representative BBH-puncture grid (the "faked" analytic initial
 * data, refined the same way the solver does) or a clean uniform grid, then runs
 * a fixed number of RK steps with the FULL per-step compute (ghost exchange,
 * unzip, derivatives, RHS, zip, constraints) but with NO remeshing and NO file
 * IO. It emits the same per-step profiling JSONL as bssnSolver so the existing
 * closure-analysis tooling works unchanged.
 *
 * It is built through the same bssn_add_executable / bssn_add_def flag system as
 * bssnSolver, so a given compile-flag combo's binary reflects exactly what the
 * solver runs (DENDRO_HYBRID_OMP, DENDRO_USE_NEW_DERIVS, cascade AVX, ...).
 *
 * Usage:
 *   bssnScalingBench <paramFile> [--mode strong|weak] [--grid bbh|uniform]
 *                    [--steps N] [--warmup K] [--lev L] [--prefix P]
 * Sweep ranks (-np) and OMP_NUM_THREADS externally; each run writes
 *   <prefix>_steps.jsonl with one record per timed step.
 */

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "TreeNode.h"
#include "bssnCtx.h"
#include "gr.h"
#include "grUtils.h"
#include "mathUtils.h"
#include "mesh.h"
#include "meshUtils.h"
#include "mpi.h"
#include "octUtils.h"
#include "parameters.h"

namespace {

struct BenchOpts {
    std::string param_file;
    std::string mode   = "strong";  // strong | weak
    std::string grid   = "bbh";     // bbh | uniform
    unsigned int steps = 10;
    unsigned int warmup = 2;
    unsigned int lev   = 3;  // uniform-grid refinement level (8^lev elements)
    // Initial field data: -1 = use the param file's BSSN_ID_TYPE (default;
    // = puncture, the only valid/stable BSSN state). Override with --id-type to
    // experiment (e.g. on a uniform grid). NOTE: a uniform grid with puncture
    // data is only marginally stable (the singularity isn't resolved) -- the
    // representative `bbh` grid is the robust multi-step benchmark.
    int id_type        = -1;
    std::string prefix = "bssn_bench";
};

void usage(const char* a0) {
    std::cout
        << "Usage: " << a0 << " <paramFile> [options]\n"
        << "  --mode strong|weak   strong = fixed grid (vary ranks); weak = fix "
           "work/rank (default strong)\n"
        << "  --grid bbh|uniform   bbh = representative puncture grid; uniform = "
           "clean kernel scaling (default bbh)\n"
        << "  --steps N            timed RK steps (default 10)\n"
        << "  --warmup K           untimed warm-up RK steps (default 2)\n"
        << "  --lev L              uniform-grid depth, 8^L elements (default 3)\n"
        << "  --id-type N          override BSSN_ID_TYPE init data (default: "
           "use param file; 0/1 = puncture)\n"
        << "  --prefix P           output prefix for <P>_steps.jsonl (default "
           "bssn_bench)\n";
}

BenchOpts parse(int argc, char** argv) {
    BenchOpts o;
    if (argc >= 2) o.param_file = argv[1];
    for (int i = 2; i < argc; i++) {
        const std::string a = argv[i];
        auto need           = [&](const char* name) -> std::string {
            if (i + 1 >= argc) {
                std::cerr << name << " requires a value\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            return argv[++i];
        };
        if (a == "--mode")
            o.mode = need("--mode");
        else if (a == "--grid")
            o.grid = need("--grid");
        else if (a == "--steps")
            o.steps = (unsigned int)std::stoul(need("--steps"));
        else if (a == "--warmup")
            o.warmup = (unsigned int)std::stoul(need("--warmup"));
        else if (a == "--lev")
            o.lev = (unsigned int)std::stoul(need("--lev"));
        else if (a == "--id-type")
            o.id_type = std::stoi(need("--id-type"));
        else if (a == "--prefix")
            o.prefix = need("--prefix");
        else
            std::cerr << "[bench] ignoring unknown arg: " << a << "\n";
    }
    return o;
}

}  // namespace

int main(int argc, char** argv) {
    // Hybrid path runs OpenMP inside compute regions, not around MPI calls.
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    const BenchOpts opt = parse(argc, argv);
    if (opt.param_file.empty()) {
        if (!rank) usage(argv[0]);
        MPI_Finalize();
        return 0;
    }

    // ---- parameters (same loader as the solver) ------------------------
    bssn::readParamFile(opt.param_file.c_str(), comm);
    if (opt.id_type >= 0) bssn::BSSN_ID_TYPE = (unsigned int)opt.id_type;
    _InitializeHcurve(bssn::BSSN_DIM);
    m_uiMaxDepth                   = bssn::BSSN_MAXDEPTH;
    bssn::BSSN_PROFILE_FILE_PREFIX = opt.prefix;

    if (!rank)
        std::cout << "[bench] mode=" << opt.mode << " grid=" << opt.grid
                  << " steps=" << opt.steps << " warmup=" << opt.warmup
                  << " lev=" << opt.lev << " npes=" << npes
                  << " rk_type=" << bssn::BSSN_RK_TYPE << " (no remesh, no IO)"
                  << std::endl;

    // ---- build the grid ------------------------------------------------
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double, double, double, double*)> f_init =
        [](double x, double y, double z, double* var) {
            bssn::punctureData(x, y, z, var);
        };
    const unsigned int interpVars = bssn::BSSN_NUM_VARS;
    unsigned int varIndex[interpVars];
    for (unsigned int i = 0; i < interpVars; i++) varIndex[i] = i;

    if (opt.grid == "uniform") {
        // Clean kernel-scaling grid: uniform octree of depth lev_use. For weak
        // scaling we grow the depth with rank count so #elements/rank stays
        // ~constant (8^lev per rank).
        unsigned int lev_use = opt.lev;
        if (opt.mode == "weak")
            lev_use += (unsigned int)std::lround(std::log((double)npes) /
                                                 std::log(8.0));
        createRegularOctree(tmpNodes, lev_use, bssn::BSSN_DIM, m_uiMaxDepth,
                            comm);
        if (!rank)
            std::cout << "[bench] uniform grid: lev=" << lev_use << std::endl;
    } else {
        // Representative BBH grid from the puncture ("faked") initial data,
        // refined exactly like the solver's adaptive init.
        const unsigned int f2olmin =
            std::min(bssn::BSSN_BH1_MAX_LEV, bssn::BSSN_BH2_MAX_LEV);
        if (f2olmin < MAXDEAPTH_LEVEL_DIFF + 2) {
            if (!rank)
                std::cout << "BH min level must be > " << (MAXDEAPTH_LEVEL_DIFF + 2)
                          << std::endl;
            MPI_Abort(comm, 0);
        }
        function2Octree(f_init, bssn::BSSN_NUM_VARS, varIndex, interpVars,
                        tmpNodes, (f2olmin - MAXDEAPTH_LEVEL_DIFF - 2),
                        bssn::BSSN_WAVELET_TOL, bssn::BSSN_ELE_ORDER, comm);
    }

    ot::Mesh* mesh = ot::createMesh(
        tmpNodes.data(), tmpNodes.size(), bssn::BSSN_ELE_ORDER, comm, 1,
        ot::SM_TYPE::FDM, bssn::BSSN_DENDRO_GRAIN_SZ, bssn::BSSN_LOAD_IMB_TOL,
        bssn::BSSN_SPLIT_FIX);
    mesh->setDomainBounds(
        Point(bssn::BSSN_GRID_MIN_X, bssn::BSSN_GRID_MIN_Y, bssn::BSSN_GRID_MIN_Z),
        Point(bssn::BSSN_GRID_MAX_X, bssn::BSSN_GRID_MAX_Y, bssn::BSSN_GRID_MAX_Z));
    unsigned int lmin, lmax;
    mesh->computeMinMaxLevel(lmin, lmax);
    tmpNodes.clear();

    // Disable AMR / remeshing: the benchmark runs on a FIXED grid (this is also
    // what keeps BSSNCtx::initialize() from doing a collective init grid-converge
    // remesh loop, which otherwise hangs/varies at multi-rank). Mirrors
    // gr_scaling.cpp. The grid we just built is what gets timed.
    bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY = 1;

    // dx / dt (same expressions as the solver)
    bssn::BSSN_CURRENT_MIN_DX =
        ((bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
         ((1u << (m_uiMaxDepth - lmax)) / ((double)bssn::BSSN_ELE_ORDER)) /
         ((double)(1u << (m_uiMaxDepth))));
    bssn::BSSN_RK45_TIME_STEP_SIZE =
        bssn::BSSN_CFL_FACTOR * bssn::BSSN_CURRENT_MIN_DX;

    // weak scaling on the BBH grid: split so each rank keeps ~grain-size elems.
    ot::Mesh* pMesh = mesh;
    if (opt.mode == "weak" && opt.grid == "bbh")
        pMesh = bssn::weakScalingReMesh(mesh, npes);

    // ---- ctx + stepper (mirror the solver, incl. MSRK) -----------------
    bssn::BSSNCtx* bssnCtx = new bssn::BSSNCtx(pMesh);
    const RKType rkType    = (RKType)bssn::BSSN_RK_TYPE;
    ts::ETS<DendroScalar, bssn::BSSNCtx>* ets = nullptr;
    if (rkType == RKType::RK4_MSRK2_1 || rkType == RKType::RK4_MSRK2_2 ||
        rkType == RKType::RK4_MSRK3) {
        ts::ETSType v = (rkType == RKType::RK4_MSRK2_1) ? ts::ETSType::RK4_MSRK2_1
                        : (rkType == RKType::RK4_MSRK2_2)
                            ? ts::ETSType::RK4_MSRK2_2
                            : ts::ETSType::RK4_MSRK3;
        ets = new ts::ETS_MSRK<DendroScalar, bssn::BSSNCtx>(bssnCtx, v);
    } else {
        ets = new ts::ETS<DendroScalar, bssn::BSSNCtx>(bssnCtx);
        if (rkType == RKType::RK3)
            ets->set_ets_coefficients(ts::ETSType::RK3);
        else if (rkType == RKType::RK5)
            ets->set_ets_coefficients(ts::ETSType::RK5);
        else
            ets->set_ets_coefficients(ts::ETSType::RK4);
    }
    ets->set_evolve_vars(bssnCtx->get_evolution_vars());
    ets->init();

    {
        // Report the actual element count on the ctx's (post-init) mesh.
        const ot::Mesh* cmesh = ets->get_mesh();
        DendroIntL le = cmesh->isActive() ? cmesh->getNumLocalMeshElements() : 0;
        DendroIntL ge = 0;
        par::Mpi_Reduce(&le, &ge, 1, MPI_SUM, 0, comm);
        if (!rank)
            std::cout << "[bench] total elements=" << ge << " lmin=" << lmin
                      << " lmax=" << lmax << std::endl;
    }

    // ---- warm-up + timed loop: NO remesh, NO IO ------------------------
    const DendroIntL start_step = ets->curr_step();
    const DendroIntL last_step  = start_step + opt.warmup + opt.steps;
    while (ets->curr_step() < last_step) {
        const DendroIntL step = ets->curr_step();
        const DendroIntL rel  = step - start_step;

        bssn::BSSN_CURRENT_RK_COORD_TIME = ets->curr_time();
        bssn::BSSN_CURRENT_RK_STEP       = step;

        // Zero every timer right before the first TIMED step (discard warm-up).
        if (rel == opt.warmup) {
#if defined(__PROFILE_ETS__) && defined(__PROFILE_CTX__)
            ets->reset_pt();
            bssnCtx->reset_pt();
#endif
            bssn::timer::resetSnapshot();
        }

        ets->evolve();

        // One JSONL record per timed step; reset after so each snap is one step.
        if (rel >= opt.warmup) {
#ifdef ENABLE_DENDRO_PROFILE_COUNTERS
#if defined(__PROFILE_ETS__) && defined(__PROFILE_CTX__)
            bssn::timer::profileInfoJSON(bssn::BSSN_PROFILE_FILE_PREFIX.c_str(),
                                         ets->get_mesh(), (unsigned int)step,
                                         &ets->m_uiCtxpt, &bssnCtx->m_uiCtxpt);
            ets->reset_pt();
#else
            bssn::timer::profileInfoJSON(bssn::BSSN_PROFILE_FILE_PREFIX.c_str(),
                                         ets->get_mesh(), (unsigned int)step);
#endif
            bssn::timer::resetSnapshot();
#endif
        }

        if (!rank)
            std::cout << "[bench] step " << step
                      << (rel < opt.warmup ? " (warmup)" : " (timed)")
                      << " t=" << ets->curr_time() << std::endl;
    }

    if (!rank)
        std::cout << "[bench] done: " << opt.steps << " timed steps -> "
                  << opt.prefix << "_steps.jsonl" << std::endl;

    // ---- cleanup -------------------------------------------------------
    ot::Mesh* final_mesh = bssnCtx->get_mesh();
    delete bssnCtx;
    delete final_mesh;
    delete ets;

    MPI_Finalize();
    return 0;
}

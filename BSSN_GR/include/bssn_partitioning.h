#pragma once

// Graph-repartitioning integration for BSSN/dendrogr.
//
// Ported from the EM4 solver's em4_partitioning.h. When
// BSSN_WITH_GRAPH_PARTITIONING is defined at compile time, these helpers
// let the solver run the ETS evolution on a graph-partitioned mesh
// (fastpart) chosen by BSSN_PARTITIONING_METHOD at runtime. When the
// compile flag is OFF, the helpers compile to no-ops and the solver
// behaves identically to the baseline SFC path.
//
// Two entry points are used:
//   - maybe_apply_partitioning() (post-IC, in bssngr_main): swap the SFC
//     initial mesh to a graph twin once the initial conditions are set.
//   - the BSSNCtx::remesh_and_gridtransfer override: the SFC<->graph
//     "sandwich" around each AMR remesh (ReMesh + interGridTransfer assume
//     SFC ordering, so we convert to an SFC twin, remesh+transfer there,
//     then convert the result back to a graph partition).
//
// A/B testing recipe:
//   - Build with BSSN_WITH_GRAPH_PARTITIONING=OFF                -> SFC
//   - Build with BSSN_WITH_GRAPH_PARTITIONING=ON, set
//     BSSN_PARTITIONING_METHOD=3 in the par file                 -> fastpart
//   - Compare per-step norms / constraints between the two.

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "dvec.h"
#include "mesh.h"
#include "meshUtils.h"
#include "parameters.h"

namespace bssn {

#ifdef BSSN_WITH_GRAPH_PARTITIONING

// In-place graph repartition of `mesh`. No-op when the runtime method is 0
// (matches SFC baseline). Returns true if a real repartition happened.
// Per-var phantom-cg fill: FLAT SPACE. Graph-partition phantom positions
// (no canonical writer / no redistribute delivery) get these instead of
// the mesh default 0 — a zero metric (det(gt)=0) NaNs the algebraic
// constraint enforcement in post_timestep (2026-06-12 BBH tracker-NaN
// root cause).
inline void registerPhantomFill(ot::Mesh* mesh) {
    std::vector<double> pf(bssn::BSSN_NUM_VARS, 0.0);
    pf[VAR::U_ALPHA]  = 1.0;
    pf[VAR::U_CHI]    = 1.0;
    pf[VAR::U_SYMGT0] = 1.0;  // gt_xx
    pf[VAR::U_SYMGT3] = 1.0;  // gt_yy
    pf[VAR::U_SYMGT5] = 1.0;  // gt_zz
    mesh->setPhantomFillValues(std::move(pf));
}

inline bool applyPartitioning(ot::Mesh* mesh) {
    if (mesh == nullptr) return false;
    if (!mesh->isActive()) return false;
    if (bssn::BSSN_PARTITIONING_METHOD == 0) return false;
    auto method =
        static_cast<PartitioningOptions>(bssn::BSSN_PARTITIONING_METHOD);
    mesh->setPartitioningMethod(method);
    mesh->repartitionMeshGlobal();
    registerPhantomFill(mesh);
    return true;
}

// Build a graph-partitioned twin of `srcMesh` (same global element set, new
// rank distribution). Caller owns the returned mesh.
//
// Injects srcMesh's ownership masks and canonical block info into the twin
// BEFORE repartitionMeshGlobal so one-sided FD stencils and block-boundary
// decomposition stay partition-invariant. Both meshes hold the same global
// element set, SFC-sorted locally, so srcMesh's local elements at
// [srcLB,srcLE) correspond positionally to the twin's [dstLB,...).
inline ot::Mesh* buildGraphTwin(ot::Mesh* srcMesh, unsigned int eleOrder,
                                unsigned int grainSz, double ldTol,
                                unsigned int splitFix, const Point& pt_min,
                                const Point& pt_max) {
    MPI_Comm comm = srcMesh->getMPIGlobalCommunicator();

    // ranks outside the active comm cannot join the labelled rebuild
    if (!srcMesh->allRanksActive()) return srcMesh;

    // DENDRO_S2G_DIRECT (default on, =0 reverts): no E2E_ONLY twin; the labelled
    // mesh is built from srcMesh's own E2E / element-ghost arrays
    // (Mesh::repartitionFrom). _CHECK=1 also runs the twin path and aborts on
    // any difference. BSSN keeps the consensus pos-bcast (repartitionMeshGlobal
    // enables it).
    static const bool s2gDirect = [] {
        const char* v = std::getenv("DENDRO_S2G_DIRECT");
        return !(v && v[0] == '0');
    }();
    static const bool s2gCheck = [] {
        const char* v = std::getenv("DENDRO_S2G_DIRECT_CHECK");
        return v && v[0] == '1';
    }();
    static int s2gCall   = 0;
    ot::Mesh* directMesh = nullptr;
    if (s2gDirect) {
        const double td0 = MPI_Wtime();
        directMesh       = ot::Mesh::repartitionFrom(
            *srcMesh,
            static_cast<PartitioningOptions>(bssn::BSSN_PARTITIONING_METHOD),
            pt_min, pt_max);
        const double td = MPI_Wtime() - td0;
        int rk          = 0;
        MPI_Comm_rank(comm, &rk);
        if (!rk)
            std::cout << "[s2g-direct] rebuild " << s2gCall << ": "
                      << (directMesh ? "DIRECT (no E2E twin)"
                                     : "fallback to twin (src not cloneable)")
                      << " " << td * 1e3 << " ms" << std::endl;
        if (directMesh && !s2gCheck) {
            s2gCall++;
            registerPhantomFill(directMesh);
            return directMesh;
        }
    }
    const bool s2gDoCheck = s2gCheck && directMesh != nullptr;
    const int s2gCallNo   = s2gCall++;
    std::vector<ot::TreeNode> oct;
    if (srcMesh->isActive()) {
        const auto* pN = srcMesh->getAllElements().data();
        for (unsigned int e = srcMesh->getElementLocalBegin();
             e < srcMesh->getElementLocalEnd(); e++)
            oct.push_back(pN[e]);
    }
    ot::Mesh* graphMesh =
        ot::createMesh(oct.data(), oct.size(), eleOrder, comm, 1,
                       ot::SM_TYPE::E2E_ONLY, grainSz, ldTol, splitFix);
    graphMesh->setDomainBounds(pt_min, pt_max);
    graphMesh->setPartitioningMethod(
        static_cast<PartitioningOptions>(bssn::BSSN_PARTITIONING_METHOD));
    graphMesh->setScatterMapType(ot::SM_TYPE::FDM);
    // consensus pos-bcast: on by default for repartitioned meshes (bssn
    // keeps the full machinery).

    if (srcMesh->isActive() && graphMesh->isActive()) {
        const unsigned int srcLB = srcMesh->getElementLocalBegin();
        const unsigned int srcLE = srcMesh->getElementLocalEnd();
        const unsigned int dstLB = graphMesh->getElementLocalBegin();
        const size_t totalElems  = graphMesh->getAllElements().size();
        const auto& srcBlkInfo = srcMesh->getBlockInfo();
        std::vector<ot::Mesh::CanonicalBlockInfo> dstBlkInfo(totalElems);
        if (srcBlkInfo.size() >= srcLE) {
            for (unsigned int i = 0;
                 i < (srcLE - srcLB) && (dstLB + i) < totalElems; i++)
                dstBlkInfo[dstLB + i] = srcBlkInfo[srcLB + i];
        }
        graphMesh->setBlockInfo(std::move(dstBlkInfo));
    }

    graphMesh->repartitionMeshGlobal();
    // _CHECK: direct labelled mesh vs the twin path's, key arrays
    if (s2gDoCheck) {
        std::string d = directMesh->diffKeyArrays(*graphMesh, false);
        int bad = d.empty() ? 0 : 1, badAll = 0;
        MPI_Allreduce(&bad, &badAll, 1, MPI_INT, MPI_SUM, comm);
        int rk = 0;
        MPI_Comm_rank(comm, &rk);
        if (bad)
            std::cerr << "[s2g-direct-check] rank " << rk << " MISMATCH in " << d
                      << std::endl;
        if (!rk)
            std::cout << "[s2g-direct-check] rebuild " << s2gCallNo << ": "
                      << badAll << " rank mismatches" << std::endl;
        if (badAll) MPI_Abort(comm, 77);
        delete directMesh;
    }
    registerPhantomFill(graphMesh);
    return graphMesh;
}

// Redistribute a DVec across a partition change. Both meshes must share the
// same global element set. On return `dvec` is bound to a fresh buffer on
// dstMesh.
//
// Fix B (from EM4, 2026-05-15): skip the src/dst performGhostExchange. The
// ETS already canonicalized every LOCAL cg at consensus phys positions at
// end-of-step; re-running performGhostExchange on a graph mesh would mirror
// stale values from paired GHOST cgs. The post-redistribute
// broadcastCgValuesByPhysPosPublic (gated by DENDRO_FORCE_POS_BCAST) brings
// all consensus positions to bit-identity instead.
template <typename T, typename I>
inline void redistributeDVec(ot::Mesh* srcMesh, ot::Mesh* dstMesh,
                             ot::DVector<T, I>& dvec) {
    ot::DVector<T, I> dvecNew;
    dvecNew.create_vector(dstMesh, dvec.get_type(), dvec.get_loc(),
                          dvec.get_dof(), dvec.is_ghost_allocated());

    if (!srcMesh->isActive() && !dstMesh->isActive()) {
        std::swap(dvec, dvecNew);
        dvecNew.destroy_vector();
        return;
    }

    const unsigned int dof = dvec.get_dof();
    T* srcPtr              = dvec.get_vec_ptr();
    T* dstPtr              = dvecNew.get_vec_ptr();
    const I sz_per_dof_old = (dof != 0) ? dvec.get_size() / dof : 0;
    const I sz_per_dof_new = (dof != 0) ? dvecNew.get_size() / dof : 0;

    // OCT_SHARED_NODES (CG) is the only redistribution path; unzip / DG /
    // cell-centered DVecs are recreated from scratch by the caller.
    if (dvec.get_type() != ot::DVEC_TYPE::OCT_SHARED_NODES) {
        std::cerr << "[bssn_partitioning] redistributeDVec: unsupported "
                     "DVec type. Skipping.\n";
        std::swap(dvec, dvecNew);
        dvecNew.destroy_vector();
        return;
    }

    // single batched redistribute: the scatter plan is built once inside and
    // all `dof` fields move together in one set of collectives (was one
    // Allgatherv+Alltoallv set per field). fix B: no ghost exchange here.
    (void)sz_per_dof_old;
    (void)sz_per_dof_new;
    srcMesh->redistributeVec(dstMesh, srcPtr, dstPtr, dof);

    // post-redistribute consensus broadcast (gated by DENDRO_FORCE_POS_BCAST
    // inside the method). env DENDRO_DISABLE_REDIST_BCAST=1 skips it (A/B).

    static const char* drb_env = std::getenv("DENDRO_DISABLE_REDIST_BCAST");
    static const bool drb_on =
        drb_env && drb_env[0] == '1' && drb_env[1] == '\0';
    if (!drb_on && dstMesh->getMPICommSize() > 1 && dstMesh->isActive()) {
        dstMesh->broadcastCgValuesByPhysPosPublic(dstPtr, dof);
    }


    std::swap(dvec, dvecNew);
    dvecNew.destroy_vector();
}

// Build an SFC-partitioned twin of `srcMesh` (same global element set, SFC
// rank distribution). Used to convert a graph mesh back to SFC before
// ReMesh / interGridTransfer (both assume SFC local-element ordering).
// Uses SM_TYPE::FDM (E2E_ONLY leaves m_uiNumActualNodes uninitialized and
// the downstream ReMesh constructor reads the source nodal scatter state).
inline ot::Mesh* buildSFCTwin(ot::Mesh* srcMesh, unsigned int eleOrder,
                              unsigned int grainSz, double ldTol,
                              unsigned int splitFix, const Point& pt_min,
                              const Point& pt_max, bool blockSetup = false) {
    MPI_Comm comm = srcMesh->getMPIGlobalCommunicator();
    std::vector<ot::TreeNode> oct;
    if (srcMesh->isActive()) {
        const auto* pN = srcMesh->getAllElements().data();
        for (unsigned int e = srcMesh->getElementLocalBegin();
             e < srcMesh->getElementLocalEnd(); e++)
            oct.push_back(pN[e]);
    }
    // blockSetup=true: a full mesh the solver can step and test for remesh on.
    // blockSetup=false: the twin is transient and never unzips — block
    // decomposition + zip plan are ~90% of the FDM ctor cost and unused.
    // (A/B'd 2026-06-11: BBH collapse reproduced with full-block twin +
    // forced twin plan, so the twin is NOT the cause.)
    ot::Mesh* sfcMesh = ot::createMesh(
        oct.data(), oct.size(), eleOrder, comm, 1, ot::SM_TYPE::FDM, grainSz,
        ldTol, splitFix, /*getWeight*/ NULL, blockSetup);
    sfcMesh->setDomainBounds(pt_min, pt_max);
    // re-arm the flag so ReMesh's successor mesh builds its blocks (the
    // twin itself still never unzips; no block state is faked).
    sfcMesh->setBlockSetupFlag(true);
    // no setPartitioningMethod / repartitionMeshGlobal -- stays SFC.
    return sfcMesh;
}

// Propagate refinement flags from srcMesh's local TreeNodes to dstMesh's
// local TreeNodes. The meshes share the global element set but differ in
// rank assignment, so allgather all local (TreeNode, flag) pairs from src,
// then each rank scans its dst local elements for matches. Applied via
// setMeshRefinementFlags so the COARSE sibling-contiguity check is enforced.
inline void redistributeFlags(ot::Mesh* srcMesh, ot::Mesh* dstMesh) {
    MPI_Comm comm = srcMesh->getMPIGlobalCommunicator();
    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    std::vector<ot::TreeNode> myLocal;
    if (srcMesh->isActive()) {
        const auto& src = srcMesh->getAllElements();
        for (unsigned int e = srcMesh->getElementLocalBegin();
             e < srcMesh->getElementLocalEnd(); e++)
            myLocal.push_back(src[e]);
    }
    int myCount = (int)myLocal.size();

    std::vector<int> counts(npes), offs(npes, 0);
    MPI_Allgather(&myCount, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
    int total = 0;
    for (int p = 0; p < npes; p++) {
        offs[p] = total;
        total += counts[p];
    }
    std::vector<ot::TreeNode> allLocal(total);
    MPI_Allgatherv(myLocal.data(), myCount,
                   par::Mpi_datatype<ot::TreeNode>::value(), allLocal.data(),
                   counts.data(), offs.data(),
                   par::Mpi_datatype<ot::TreeNode>::value(), comm);

    std::map<ot::TreeNode, unsigned int> tnFlag;
    for (const auto& tn : allLocal) tnFlag[tn] = tn.getFlag();

    if (!dstMesh->isActive()) return;
    const unsigned int numLocal =
        dstMesh->getElementLocalEnd() - dstMesh->getElementLocalBegin();
    std::vector<unsigned int> refine_flags(numLocal, OCT_NO_CHANGE);
    for (unsigned int e = dstMesh->getElementLocalBegin();
         e < dstMesh->getElementLocalEnd(); e++) {
        const ot::TreeNode& dst = dstMesh->getAllElements()[e];
        auto it                 = tnFlag.find(dst);
        if (it == tnFlag.end()) continue;
        const unsigned int f   = (it->second) >> NUM_LEVEL_BITS;
        const unsigned int rid = e - dstMesh->getElementLocalBegin();
        if (f == OCT_SPLIT)
            refine_flags[rid] = OCT_SPLIT;
        else if (f == OCT_COARSE)
            refine_flags[rid] = OCT_COARSE;
        // else leave OCT_NO_CHANGE
    }
    dstMesh->setMeshRefinementFlags(refine_flags);
}

#else  // !BSSN_WITH_GRAPH_PARTITIONING

inline bool applyPartitioning(ot::Mesh* /*mesh*/) { return false; }

inline ot::Mesh* buildGraphTwin(ot::Mesh* srcMesh, unsigned int /*eleOrder*/,
                                unsigned int /*grainSz*/, double /*ldTol*/,
                                unsigned int /*splitFix*/,
                                const Point& /*pt_min*/,
                                const Point& /*pt_max*/) {
    return srcMesh;
}

inline ot::Mesh* buildSFCTwin(ot::Mesh* srcMesh, unsigned int /*eleOrder*/,
                              unsigned int /*grainSz*/, double /*ldTol*/,
                              unsigned int /*splitFix*/,
                              const Point& /*pt_min*/,
                              const Point& /*pt_max*/,
                              bool /*blockSetup*/ = false) {
    return srcMesh;
}

inline void redistributeFlags(ot::Mesh* /*srcMesh*/, ot::Mesh* /*dstMesh*/) {}

template <typename T, typename I>
inline void redistributeDVec(ot::Mesh* /*srcMesh*/, ot::Mesh* /*dstMesh*/,
                             ot::DVector<T, I>& /*dvec*/) {}

#endif  // BSSN_WITH_GRAPH_PARTITIONING

}  // namespace bssn

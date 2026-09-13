#include "OnMeshTeukolskyID.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include "derivs.h"
#include "grDef.h"
#include "parameters.h"
#ifdef BUILD_WITH_PETSC
#include <petscksp.h>
#endif

namespace bssn {
namespace {

void first(unsigned axis, double* out, const double* in, double h,
           const unsigned* sz, unsigned flag) {
#ifdef DENDRO_USE_NEW_DERIVS
    if (axis == 0) BSSN_DERIVS->deriv_x(out, in, h, sz, flag);
    if (axis == 1) BSSN_DERIVS->deriv_y(out, in, h, sz, flag);
    if (axis == 2) BSSN_DERIVS->deriv_z(out, in, h, sz, flag);
#else
    const auto f = axis == 0 ? deriv_x : axis == 1 ? deriv_y : deriv_z;
    f(out, in, h, sz, flag);
#endif
}
void second(unsigned axis, double* out, const double* in, double h,
            const unsigned* sz, unsigned flag) {
#ifdef DENDRO_USE_NEW_DERIVS
    if (axis == 0) BSSN_DERIVS->deriv_xx(out, in, h, sz, flag);
    if (axis == 1) BSSN_DERIVS->deriv_yy(out, in, h, sz, flag);
    if (axis == 2) BSSN_DERIVS->deriv_zz(out, in, h, sz, flag);
#else
    const auto f = axis == 0 ? deriv_xx : axis == 1 ? deriv_yy : deriv_zz;
    f(out, in, h, sz, flag);
#endif
}

struct Norms {
    double l2 = 0.0, maximum = 0.0;
};

Norms norms(const ot::Mesh& mesh, const double* v) {
    double sum = 0.0, maximum = 0.0;
    unsigned long long count = 0;
    for (unsigned int n = mesh.getNodeLocalBegin(); n < mesh.getNodeLocalEnd();
         ++n) {
        sum += v[n] * v[n];
        maximum = std::max(maximum, std::abs(v[n]));
        ++count;
    }
    double global_sum = 0.0, global_max = 0.0;
    unsigned long long global_count = 0;
    MPI_Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM,
                  mesh.getMPICommunicator());
    MPI_Allreduce(&maximum, &global_max, 1, MPI_DOUBLE, MPI_MAX,
                  mesh.getMPICommunicator());
    MPI_Allreduce(&count, &global_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                  mesh.getMPICommunicator());
    return {
        std::sqrt(global_sum / std::max<unsigned long long>(global_count, 1)),
        global_max};
}

double dot(const ot::Mesh& mesh, const double* a, const double* b) {
    double local = 0.0, global = 0.0;
    for (unsigned int n = mesh.getNodeLocalBegin(); n < mesh.getNodeLocalEnd();
         ++n)
        local += a[n] * b[n];
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM,
                  mesh.getMPICommunicator());
    return global;
}

class Operator {
   public:
    Operator(ot::Mesh& mesh, double** vars, const double* ricci) : mesh_(mesh) {
        for (auto& p : coef_) p = mesh_.createVector<double>(0.0);
        for (auto& p : coef_uz_) p = mesh_.createUnZippedVector<double>(0.0);
        u_uz_                    = mesh_.createUnZippedVector<double>(0.0);
        out_uz_                  = mesh_.createUnZippedVector<double>(0.0);
        const unsigned int begin = mesh_.getNodeLocalBegin(),
                           end   = mesh_.getNodeLocalEnd();
        for (unsigned int n = begin; n < end; ++n) {
            const double a = vars[VAR::U_SYMGT0][n], b = vars[VAR::U_SYMGT1][n];
            const double c = vars[VAR::U_SYMGT2][n], d = vars[VAR::U_SYMGT3][n];
            const double e = vars[VAR::U_SYMGT4][n], f = vars[VAR::U_SYMGT5][n];
            const double det =
                a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d);
            coef_[0][n] = (d * f - e * e) / det;
            coef_[1][n] = (c * e - b * f) / det;
            coef_[2][n] = (b * e - c * d) / det;
            coef_[3][n] = (a * f - c * c) / det;
            coef_[4][n] = (b * c - a * e) / det;
            coef_[5][n] = (a * d - b * b) / det;
            coef_[6][n] = ricci[n];
            for (unsigned a = 0; a < 3; ++a)
                coef_[7 + a][n] = vars[VAR::U_GT0 + a][n];
        }
        for (unsigned int s = 0; s < 10; ++s) {
            mesh_.readFromGhostBegin(coef_[s], 1);
            mesh_.readFromGhostEnd(coef_[s], 1);
            mesh_.unzip(coef_[s], coef_uz_[s]);
        }
        find_boundary();
    }
    ~Operator() {
        for (auto* p : coef_) delete[] p;
        for (auto* p : coef_uz_) delete[] p;
        delete[] u_uz_;
        delete[] out_uz_;
    }
    const std::vector<unsigned int>& boundary() const { return boundary_; }
    void apply(const double* u, double* out) {
        mesh_.readFromGhostBegin(const_cast<double*>(u), 1);
        mesh_.readFromGhostEnd(const_cast<double*>(u), 1);
        mesh_.unzip(u, u_uz_);
        std::fill(out_uz_, out_uz_ + mesh_.getDegOfFreedomUnZip(), 0.0);
        const Point lo = mesh_.getDomainMinPt(), hi = mesh_.getDomainMaxPt();
        for (const auto& block : mesh_.getLocalBlockList()) {
            const unsigned int off   = block.getOffset();
            const unsigned int sz[3] = {block.getAllocationSzX(),
                                        block.getAllocationSzY(),
                                        block.getAllocationSzZ()};
            const unsigned int n     = sz[0] * sz[1] * sz[2],
                               flag  = block.getBlkNodeFlag();
            const double h[3]        = {block.computeDx(lo, hi),
                                        block.computeDy(lo, hi),
                                        block.computeDz(lo, hi)};
            std::array<std::vector<double>, 3> du, d2;
            std::array<std::vector<double>, 3> mix;
            for (auto& v : du) v.resize(n);
            for (auto& v : d2) v.resize(n);
            for (auto& v : mix) v.resize(n);
            first(0, du[0].data(), u_uz_ + off, h[0], sz, flag);
            first(1, du[1].data(), u_uz_ + off, h[1], sz, flag);
            first(2, du[2].data(), u_uz_ + off, h[2], sz, flag);
            second(0, d2[0].data(), u_uz_ + off, h[0], sz, flag);
            second(1, d2[1].data(), u_uz_ + off, h[1], sz, flag);
            second(2, d2[2].data(), u_uz_ + off, h[2], sz, flag);
            // Mixed derivatives use the configured production first-derivative
            // operator successively, matching the generated BSSN Ricci path.
            std::vector<double> scratch(n);
            first(0, scratch.data(), u_uz_ + off, h[0], sz, flag);
            first(1, mix[0].data(), scratch.data(), h[1], sz, flag);
            first(2, mix[1].data(), scratch.data(), h[2], sz, flag);
            first(1, scratch.data(), u_uz_ + off, h[1], sz, flag);
            first(2, mix[2].data(), scratch.data(), h[2], sz, flag);
            const unsigned int pw = BSSN_PADDING_WIDTH;
            for (unsigned int k = pw; k < sz[2] - pw; ++k)
                for (unsigned int j = pw; j < sz[1] - pw; ++j)
                    for (unsigned int i = pw; i < sz[0] - pw; ++i) {
                        const unsigned int p = i + sz[0] * (j + sz[1] * k);
                        const double axx     = coef_uz_[0][off + p],
                                     axy     = coef_uz_[1][off + p],
                                     axz     = coef_uz_[2][off + p];
                        const double ayy     = coef_uz_[3][off + p],
                                     ayz     = coef_uz_[4][off + p],
                                     azz     = coef_uz_[5][off + p];
                        double value         = axx * d2[0][p] + ayy * d2[1][p] +
                                       azz * d2[2][p] + 2 * axy * mix[0][p] +
                                       2 * axz * mix[1][p] +
                                       2 * ayz * mix[2][p];
                        for (unsigned a = 0; a < 3; ++a)
                            value -= coef_uz_[7 + a][off + p] * du[a][p];
                        out_uz_[off + p] = value - 0.125 *
                                                       coef_uz_[6][off + p] *
                                                       u_uz_[off + p];
                    }
        }
        mesh_.zip(out_uz_, out);
        for (auto n : boundary_) out[n] = u[n];
    }

   private:
    ot::Mesh& mesh_;
    std::array<double*, 10> coef_{}, coef_uz_{};
    double *u_uz_ = nullptr, *out_uz_ = nullptr;
    std::vector<unsigned int> boundary_;
    void find_boundary() {
        const auto& els      = mesh_.getAllElements();
        const auto& e2n      = mesh_.getE2NMapping();
        const unsigned int o = mesh_.getElementOrder(), n1 = o + 1,
                           n2 = n1 * n1, npe = n1 * n2,
                           gmax = 1u << m_uiMaxDepth;
        for (unsigned int q = mesh_.getElementLocalBegin();
             q < mesh_.getElementLocalEnd(); ++q) {
            const auto& e = els[q];
            const bool xm = e.minX() == 0, xp = e.maxX() == gmax,
                       ym = e.minY() == 0, yp = e.maxY() == gmax,
                       zm = e.minZ() == 0, zp = e.maxZ() == gmax;
            if (!(xm || xp || ym || yp || zm || zp)) continue;
            for (unsigned int k = 0; k < n1; ++k)
                for (unsigned int j = 0; j < n1; ++j)
                    for (unsigned int i = 0; i < n1; ++i)
                        if ((xm && i == 0) || (xp && i == o) ||
                            (ym && j == 0) || (yp && j == o) ||
                            (zm && k == 0) || (zp && k == o))
                            boundary_.push_back(
                                e2n[q * npe + k * n2 + j * n1 + i]);
        }
        std::sort(boundary_.begin(), boundary_.end());
        boundary_.erase(std::unique(boundary_.begin(), boundary_.end()),
                        boundary_.end());
    }
};

#ifdef BUILD_WITH_PETSC
// PETSc vectors hold only owned nodes. The shell expands them into Dendro's
// ghosted layout and delegates all mesh communication to Operator::apply.
struct PetscOperator {
    Operator& op;
    unsigned begin, end;
    std::vector<double> input, output;
};
PetscErrorCode petscMultiply(Mat matrix, Vec x, Vec y) {
    PetscOperator* ctx = nullptr;
    PetscErrorCode err = MatShellGetContext(matrix, &ctx);
    if (err) return err;
    const PetscScalar* in = nullptr;
    PetscScalar* out      = nullptr;
    err                   = VecGetArrayRead(x, &in);
    if (err) return err;
    for (unsigned n = ctx->begin; n < ctx->end; ++n)
        ctx->input[n] = PetscRealPart(in[n - ctx->begin]);
    err = VecRestoreArrayRead(x, &in);
    if (err) return err;
    ctx->op.apply(ctx->input.data(), ctx->output.data());
    err = VecGetArray(y, &out);
    if (err) return err;
    for (unsigned n = ctx->begin; n < ctx->end; ++n)
        out[n - ctx->begin] = ctx->output[n];
    return VecRestoreArray(y, &out);
}
void petscCheck(PetscErrorCode err) {
    if (err)
        throw std::runtime_error("PETSc failure in type-14 Hamiltonian solve");
}
#endif

void zero_boundary(const std::vector<unsigned int>& boundary, double* v) {
    for (auto n : boundary) v[n] = 0.0;
}

}  // namespace

std::pair<double, double> teukolskyConnectionOnMesh(ot::Mesh& mesh,
                                                    double** vars,
                                                    bool initialize) {
    if (!mesh.isActive()) return {0.0, 0.0};
    std::array<std::vector<double>, 6> metric;
    std::array<std::vector<double>, 3> output, zipped;
    const unsigned nu = mesh.getDegOfFreedomUnZip(),
                   nz = mesh.getDegOfFreedom();
    for (unsigned s = 0; s < 6; ++s) {
        metric[s].resize(nu);
        double* v = vars[VAR::U_SYMGT0 + s];
        mesh.readFromGhostBegin(v, 1);
        mesh.readFromGhostEnd(v, 1);
        mesh.unzip(v, metric[s].data());
    }
    for (unsigned a = 0; a < 3; ++a) {
        output[a].resize(nu);
        zipped[a].resize(nz);
    }
    const Point lo = mesh.getDomainMinPt(), hi = mesh.getDomainMaxPt();
    for (const auto& block : mesh.getLocalBlockList()) {
        const unsigned off   = block.getOffset();
        const unsigned sz[3] = {block.getAllocationSzX(),
                                block.getAllocationSzY(),
                                block.getAllocationSzZ()};
        const unsigned n = sz[0] * sz[1] * sz[2], flag = block.getBlkNodeFlag();
        const double h[3] = {block.computeDx(lo, hi), block.computeDy(lo, hi),
                             block.computeDz(lo, hi)};
        std::array<std::array<std::vector<double>, 3>, 6> grad;
        for (unsigned s = 0; s < 6; ++s)
            for (unsigned a = 0; a < 3; ++a) {
                grad[s][a].resize(n);
                first(a, grad[s][a].data(), metric[s].data() + off, h[a], sz,
                      flag);
            }
        const double* gt0        = metric[0].data() + off;
        const double* grad_0_gt0 = grad[0][0].data();
        const double* grad_1_gt0 = grad[0][1].data();
        const double* grad_2_gt0 = grad[0][2].data();
        const double* gt1        = metric[1].data() + off;
        const double* grad_0_gt1 = grad[1][0].data();
        const double* grad_1_gt1 = grad[1][1].data();
        const double* grad_2_gt1 = grad[1][2].data();
        const double* gt2        = metric[2].data() + off;
        const double* grad_0_gt2 = grad[2][0].data();
        const double* grad_1_gt2 = grad[2][1].data();
        const double* grad_2_gt2 = grad[2][2].data();
        const double* gt3        = metric[3].data() + off;
        const double* grad_0_gt3 = grad[3][0].data();
        const double* grad_1_gt3 = grad[3][1].data();
        const double* grad_2_gt3 = grad[3][2].data();
        const double* gt4        = metric[4].data() + off;
        const double* grad_0_gt4 = grad[4][0].data();
        const double* grad_1_gt4 = grad[4][1].data();
        const double* grad_2_gt4 = grad[4][2].data();
        const double* gt5        = metric[5].data() + off;
        const double* grad_0_gt5 = grad[5][0].data();
        const double* grad_1_gt5 = grad[5][1].data();
        const double* grad_2_gt5 = grad[5][2].data();
        double* connection0      = output[0].data() + off;
        double* connection1      = output[1].data() + off;
        double* connection2      = output[2].data() + off;
        const unsigned pw        = BSSN_PADDING_WIDTH;
        for (unsigned k = pw; k < sz[2] - pw; ++k)
            for (unsigned j = pw; j < sz[1] - pw; ++j)
                for (unsigned i = pw; i < sz[0] - pw; ++i) {
                    const unsigned pp = i + sz[0] * (j + sz[1] * k);
#include "teuk_connection.inc"
                }
    }
    for (unsigned a = 0; a < 3; ++a)
        mesh.zip(output[a].data(), zipped[a].data());
    std::vector<double> error(nz, 0.0);
    for (unsigned n = mesh.getNodeLocalBegin(); n < mesh.getNodeLocalEnd();
         ++n) {
        for (unsigned a = 0; a < 3; ++a) {
            if (initialize) vars[VAR::U_GT0 + a][n] = zipped[a][n];
            const double d = vars[VAR::U_GT0 + a][n] - zipped[a][n];
            error[n] += d * d;
        }
        error[n] = std::sqrt(error[n]);
    }
    const auto norm = norms(mesh, error.data());
    return {norm.l2, norm.maximum};
}

TeukHamiltonianSolveResult solveTeukolskyHamiltonianOnMesh(
    ot::Mesh& mesh, double** vars, const double* ricci, double tolerance,
    unsigned int maximum_iterations, bool verbose) {
    Operator op(mesh, vars, ricci);
    const unsigned int size  = mesh.getDegOfFreedom();
    const unsigned int begin = mesh.getNodeLocalBegin(),
                       end   = mesh.getNodeLocalEnd();
    auto make                = [&]() { return mesh.createVector<double>(0.0); };
    double *u = make(), *rhs = make(), *ax = make(), *r = make(),
           *shadow = make(), *p = make(), *v = make(), *s = make(), *t = make();
    for (unsigned int n = begin; n < end; ++n) rhs[n] = 0.125 * ricci[n];
    zero_boundary(op.boundary(), rhs);
    op.apply(u, ax);
    for (unsigned int n = begin; n < end; ++n) r[n] = rhs[n] - ax[n];
    zero_boundary(op.boundary(), r);
    std::copy(r, r + size, shadow);
    const double base = std::sqrt(dot(mesh, r, r));
    TeukHamiltonianSolveResult result;
    result.converged = (base == 0.0);
    double rho_old = 1, alpha = 1, omega = 1;
#ifdef BUILD_WITH_PETSC
    PetscOperator shell{op, begin, end, std::vector<double>(size),
                        std::vector<double>(size)};
    Mat matrix = nullptr;
    Vec x = nullptr, b = nullptr;
    KSP ksp             = nullptr;
    const MPI_Comm comm = mesh.getMPICommunicator();
    petscCheck(MatCreateShell(comm, end - begin, end - begin, PETSC_DETERMINE,
                              PETSC_DETERMINE, &shell, &matrix));
    petscCheck(MatShellSetOperation(
        matrix, MATOP_MULT, reinterpret_cast<void (*)(void)>(petscMultiply)));
    petscCheck(VecCreateMPI(comm, end - begin, PETSC_DETERMINE, &x));
    petscCheck(VecDuplicate(x, &b));
    petscCheck(VecSet(x, 0.0));
    PetscScalar* array = nullptr;
    petscCheck(VecGetArray(b, &array));
    for (unsigned n = begin; n < end; ++n) array[n - begin] = rhs[n];
    petscCheck(VecRestoreArray(b, &array));
    petscCheck(KSPCreate(comm, &ksp));
    petscCheck(KSPSetOperators(ksp, matrix, matrix));
    petscCheck(KSPSetType(ksp, KSPGMRES));
    petscCheck(KSPGMRESSetRestart(ksp, 100));
    PC pc = nullptr;
    petscCheck(KSPGetPC(ksp, &pc));
    petscCheck(PCSetType(pc, PCNONE));
    petscCheck(KSPSetTolerances(ksp, tolerance, 1e-50, PETSC_DEFAULT,
                                maximum_iterations));
    petscCheck(KSPSetInitialGuessNonzero(ksp, PETSC_FALSE));
    if (verbose)
        petscCheck(KSPMonitorSet(
            ksp,
            [](KSP solver, PetscInt it, PetscReal residual,
               void*) -> PetscErrorCode {
                int rank = 0;
                MPI_Comm_rank(PetscObjectComm((PetscObject)solver), &rank);
                if (!rank && (it == 0 || it % 25 == 0))
                    std::cout << "TEUK_HAM PETSc iteration " << it
                              << " residual " << residual << '\n';
                return 0;
            },
            nullptr, nullptr));
    petscCheck(KSPSolve(ksp, b, x));
    PetscInt iterations = 0;
    petscCheck(KSPGetIterationNumber(ksp, &iterations));
    result.iterations           = iterations;
    const PetscScalar* solution = nullptr;
    petscCheck(VecGetArrayRead(x, &solution));
    for (unsigned n = begin; n < end; ++n)
        u[n] = PetscRealPart(solution[n - begin]);
    petscCheck(VecRestoreArrayRead(x, &solution));
    petscCheck(KSPDestroy(&ksp));
    petscCheck(VecDestroy(&x));
    petscCheck(VecDestroy(&b));
    petscCheck(MatDestroy(&matrix));
#else
    for (unsigned int it = 1; it <= maximum_iterations && base > 0; ++it) {
        const double rho = dot(mesh, shadow, r);
        if (!std::isfinite(rho) || std::abs(rho) < 1e-300) break;
        const double beta = (rho / rho_old) * (alpha / omega);
        for (unsigned int n = begin; n < end; ++n)
            p[n] = r[n] + beta * (p[n] - omega * v[n]);
        zero_boundary(op.boundary(), p);
        op.apply(p, v);
        const double den = dot(mesh, shadow, v);
        if (!std::isfinite(den) || std::abs(den) < 1e-300) break;
        alpha = rho / den;
        for (unsigned int n = begin; n < end; ++n) s[n] = r[n] - alpha * v[n];
        zero_boundary(op.boundary(), s);
        if (std::sqrt(dot(mesh, s, s)) <= tolerance * base) {
            for (unsigned int n = begin; n < end; ++n) u[n] += alpha * p[n];
            result.iterations = it;
            result.converged  = true;
            break;
        }
        op.apply(s, t);
        const double tt = dot(mesh, t, t);
        if (!std::isfinite(tt) || tt == 0) break;
        omega = dot(mesh, t, s) / tt;
        for (unsigned int n = begin; n < end; ++n) {
            u[n] += alpha * p[n] + omega * s[n];
            r[n] = s[n] - omega * t[n];
        }
        zero_boundary(op.boundary(), r);
        result.iterations = it;
        if (std::sqrt(dot(mesh, r, r)) <= tolerance * base) {
            result.converged = true;
            break;
        }
        if (!std::isfinite(omega) || std::abs(omega) < 1e-300) break;
        rho_old = rho;
        if (verbose && (it == 1 || it % 25 == 0)) {
            const double relative = std::sqrt(dot(mesh, r, r)) / base;
            if (!mesh.getMPIRank())
                std::cout << "TEUK_HAM iteration " << it
                          << " relative residual " << relative << '\n';
        }
    }
#endif
    op.apply(u, ax);
    for (unsigned int n = begin; n < end; ++n) r[n] = ax[n] - rhs[n];
    const Norms rn = norms(mesh, r);
    result.converged =
        std::isfinite(rn.l2) && std::sqrt(dot(mesh, r, r)) <= tolerance * base;
    result.residual_l2  = rn.l2;
    result.residual_max = rn.maximum;
    double local_min = 1e300, local_max = -1e300;
    for (unsigned int n = begin; n < end; ++n) {
        const double psi = 1.0 + u[n];
        local_min =
            std::min(local_min, std::isfinite(psi)
                                    ? psi
                                    : -std::numeric_limits<double>::infinity());
        local_max = std::max(
            local_max,
            std::isfinite(psi) ? psi : std::numeric_limits<double>::infinity());
    }
    MPI_Allreduce(&local_min, &result.psi_min, 1, MPI_DOUBLE, MPI_MIN,
                  mesh.getMPICommunicator());
    MPI_Allreduce(&local_max, &result.psi_max, 1, MPI_DOUBLE, MPI_MAX,
                  mesh.getMPICommunicator());
    if (!(result.psi_min > 0.0) || !std::isfinite(result.psi_max))
        result.converged = false;
    if (result.converged)
        for (unsigned n = begin; n < end; ++n)
            vars[VAR::U_CHI][n] = std::pow(1.0 + u[n], -4.0);
    for (double* q : {u, rhs, ax, r, shadow, p, v, s, t}) delete[] q;
    return result;
}

}  // namespace bssn

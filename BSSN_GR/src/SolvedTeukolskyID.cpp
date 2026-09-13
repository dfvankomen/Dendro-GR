#include "SolvedTeukolskyID.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

#include "grDef.h"
#include "parameters.h"

namespace bssn {
namespace {

constexpr unsigned int FILE_VERSION = 1;
constexpr unsigned int FILE_FIELDS = 24;

/** Binary layout (also documented by save_dendro_binary in solve_hamiltonian.py):
 * all values are little-endian; magic[8]="DTKID13\0", uint32 version,
 * uint32 field_count, uint64 nx,ny,nz, then nine float64 values containing
 * xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,dz.  The payload is 24 field-major
 * float64 arrays in grDef.h VAR order, each indexed [ix][iy][iz] with z
 * fastest.  U_SYMGT0..5 and U_SYMAT0..5 are ordered
 * (xx,xy,xz,yy,yz,zz), as established by the type-9 assignments.
 */
class SolvedGrid {
   public:
    void evaluate(double x, double y, double z, double* var) const {
        Bracket bx, by, bz;
        if (!bracket(x, bounds_[0], spacing_[0], n_[0], bx) ||
            !bracket(y, bounds_[2], spacing_[1], n_[1], by) ||
            !bracket(z, bounds_[4], spacing_[2], n_[2], bz)) {
            asymptotic(var);
            bool expected = false;
            if (outside_warning_.compare_exchange_strong(expected, true)) {
                int rank = 0;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                if (!rank)
                    std::cerr << "WARNING: BSSN_ID_TYPE=13 requested a point outside "
                              << TEUK_SOLVED_ID_FILE
                              << "; using asymptotically flat data. First point: ("
                              << x << ", " << y << ", " << z << ")\n";
            }
            return;
        }
        for (unsigned int field = 0; field < FILE_FIELDS; ++field)
            var[field] = interpolate(field, bx, by, bz);
    }

    void load(const std::string& path) {
        const unsigned short endian_test = 1;
        if (*reinterpret_cast<const unsigned char*>(&endian_test) != 1)
            throw std::runtime_error("type-13 reader requires a little-endian host");
        if (sizeof(double) != 8 || !std::numeric_limits<double>::is_iec559)
            throw std::runtime_error("type-13 reader requires IEEE-754 float64");
        std::ifstream in(path, std::ios::binary);
        if (!in) throw std::runtime_error("cannot open solved Teukolsky file: " + path);
        char magic[8];
        read(in, magic, sizeof(magic));
        const char expected[8] = {'D','T','K','I','D','1','3','\0'};
        if (std::memcmp(magic, expected, 8) != 0)
            throw std::runtime_error("invalid type-13 magic in " + path);
        std::uint32_t version = 0, field_count = 0;
        read(in, &version, sizeof(version)); read(in, &field_count, sizeof(field_count));
        if (version != FILE_VERSION || field_count != FILE_FIELDS)
            throw std::runtime_error("unsupported type-13 version/field count in " + path);
        read(in, n_.data(), 3*sizeof(n_[0]));
        read(in, bounds_.data(), 6*sizeof(bounds_[0]));
        read(in, spacing_.data(), 3*sizeof(spacing_[0]));
        for (unsigned int d = 0; d < 3; ++d) {
            if (n_[d] < 2 || !(spacing_[d] > 0.0))
                throw std::runtime_error("invalid type-13 grid dimensions or spacing");
            const double expected_max = bounds_[2*d] + spacing_[d]*(n_[d]-1);
            const double scale = std::max(1.0, std::abs(bounds_[2*d+1]));
            if (std::abs(expected_max-bounds_[2*d+1]) > 128*std::numeric_limits<double>::epsilon()*scale)
                throw std::runtime_error("inconsistent type-13 bounds and spacing");
        }
        const size_t count = checked_count();
        for (auto& field : fields_) {
            field.resize(count);
            read(in, field.data(), count*sizeof(double));
        }
        char extra;
        if (in.read(&extra, 1))
            throw std::runtime_error("unexpected trailing bytes in type-13 file");
        int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (!rank)
            std::cout << "Loaded BSSN_ID_TYPE=13 file " << path << " grid="
                      << n_[0] << "x" << n_[1] << "x" << n_[2] << " bounds=["
                      << bounds_[0] << "," << bounds_[1] << "]x[" << bounds_[2]
                      << "," << bounds_[3] << "]x[" << bounds_[4] << ","
                      << bounds_[5] << "]\n";
        if (TEUK_SOLVED_ID_VERIFY) verify_nodes(rank);
    }

   private:
    struct Bracket { size_t lo, hi; double weight; };
    std::array<std::uint64_t, 3> n_{{0,0,0}};
    std::array<double, 6> bounds_{{0,0,0,0,0,0}};
    std::array<double, 3> spacing_{{0,0,0}};
    std::array<std::vector<double>, FILE_FIELDS> fields_;
    mutable std::atomic<bool> outside_warning_{false};

    static void read(std::ifstream& in, void* data, size_t bytes) {
        if (!in.read(reinterpret_cast<char*>(data), bytes))
            throw std::runtime_error("truncated type-13 solved Teukolsky file");
    }
    size_t checked_count() const {
        const std::uint64_t count = n_[0]*n_[1]*n_[2];
        if (n_[0] && count/n_[0]/n_[1] != n_[2] ||
            count > std::numeric_limits<size_t>::max())
            throw std::runtime_error("type-13 grid is too large");
        return static_cast<size_t>(count);
    }
    size_t index(size_t i, size_t j, size_t k) const {
        return (i*static_cast<size_t>(n_[1]) + j)*static_cast<size_t>(n_[2]) + k;
    }
    static bool bracket(double x, double minimum, double h,
                        std::uint64_t count, Bracket& b) {
        const double t = (x-minimum)/h;
        const double last = static_cast<double>(count-1);
        const double eps = 64*std::numeric_limits<double>::epsilon()*
                           std::max(1.0, std::abs(t));
        if (t < -eps || t > last+eps) return false;
        const double nearest = std::round(t);
        if (std::abs(t-nearest) <= eps) {
            b.lo = b.hi = static_cast<size_t>(std::max(0.0, std::min(last, nearest)));
            b.weight = 0.0;
            return true;
        }
        b.lo = static_cast<size_t>(std::floor(t));
        b.hi = b.lo+1;
        b.weight = t-static_cast<double>(b.lo);
        return true;
    }
    double interpolate(unsigned int f, const Bracket& x, const Bracket& y,
                       const Bracket& z) const {
        double result = 0.0;
        for (unsigned int a = 0; a < 2; ++a)
            for (unsigned int b = 0; b < 2; ++b)
                for (unsigned int c = 0; c < 2; ++c) {
                    const size_t i = a ? x.hi : x.lo, j = b ? y.hi : y.lo;
                    const size_t k = c ? z.hi : z.lo;
                    const double wx = a ? x.weight : 1.0-x.weight;
                    const double wy = b ? y.weight : 1.0-y.weight;
                    const double wz = c ? z.weight : 1.0-z.weight;
                    result += wx*wy*wz*fields_[f][index(i,j,k)];
                }
        return result;
    }
    static void asymptotic(double* var) {
        std::fill(var, var+FILE_FIELDS, 0.0);
        var[VAR::U_ALPHA] = 1.0; var[VAR::U_CHI] = 1.0;
        var[VAR::U_SYMGT0] = 1.0;  // xx
        var[VAR::U_SYMGT3] = 1.0;  // yy
        var[VAR::U_SYMGT5] = 1.0;  // zz
    }
    void verify_nodes(int rank) const {
        std::array<double, FILE_FIELDS> maximum{};
        const std::array<std::uint64_t, 3> sample_x{{0,n_[0]/2,n_[0]-1}};
        const std::array<std::uint64_t, 3> sample_y{{0,n_[1]/2,n_[1]-1}};
        const std::array<std::uint64_t, 3> sample_z{{0,n_[2]/2,n_[2]-1}};
        for (auto i : sample_x) for (auto j : sample_y) for (auto k : sample_z) {
            Bracket bx{i,i,0}, by{j,j,0}, bz{k,k,0};
            for (unsigned int f=0; f<FILE_FIELDS; ++f)
                maximum[f] = std::max(maximum[f],
                    std::abs(interpolate(f,bx,by,bz)-fields_[f][index(i,j,k)]));
        }
        if (!rank) {
            std::cout << std::setprecision(17) << "TYPE13 exact-node verification:\n"
                      << "  chi " << maximum[VAR::U_CHI] << '\n';
            for (unsigned int s=0;s<6;++s)
                std::cout << "  SYMGT" << s << " " << maximum[VAR::U_SYMGT0+s] << '\n';
            for (unsigned int s=0;s<3;++s)
                std::cout << "  GT" << s << " " << maximum[VAR::U_GT0+s] << '\n';
            for (unsigned int s=0;s<6;++s)
                std::cout << "  SYMAT" << s << " " << maximum[VAR::U_SYMAT0+s] << '\n';
            std::cout << "  K " << maximum[VAR::U_K] << std::endl;
        }
    }
};

SolvedGrid& solved_grid() {
    static SolvedGrid grid;
    static std::once_flag loaded;
    std::call_once(loaded, [&]() { grid.load(TEUK_SOLVED_ID_FILE); });
    return grid;
}

}  // namespace

void solvedTeukolskyData(const double x_grid, const double y_grid,
                         const double z_grid, double* var) {
    solved_grid().evaluate(GRIDX_TO_X(x_grid), GRIDY_TO_Y(y_grid),
                           GRIDZ_TO_Z(z_grid), var);
}

}  // namespace bssn

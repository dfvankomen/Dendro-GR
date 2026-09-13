#ifndef DENDRO_GR_SOLVED_TEUKOLSKY_ID_H
#define DENDRO_GR_SOLVED_TEUKOLSKY_ID_H

namespace bssn {

/** Interpolate file-backed, Hamiltonian-solved Teukolsky data (ID type 13).
 * Input coordinates are Dendro octree/grid coordinates, as for the other
 * initialDataFunctionWrapper targets; conversion to physical coordinates is
 * performed internally before interpolation.
 */
void solvedTeukolskyData(const double x_grid, const double y_grid,
                         const double z_grid, double* var);

}  // namespace bssn

#endif

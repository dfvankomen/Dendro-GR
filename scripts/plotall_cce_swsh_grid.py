##
# @brief: Plots the SWSH collocation grid for EVERY field in a CCE worldtube
#          HDF5 file, saving one PNG per field into an output folder. Same
#          plot (3D sphere + theta-phi map) as plot_cce_swsh_grid.py, just
#          looped over every dataset instead of picking one with --field.
#
# Usage:
#   python plotall_cce_swsh_grid.py <h5file> [--row -1] [--l-max 20]
#                                    [--radius 100] [--outdir swsh_plots]
#
# @date: 2026-10-08

import argparse
import os

import h5py
import matplotlib
matplotlib.use("Agg")  # no display needed; also avoids piling up open GUI
                        # windows across 49 figures
import matplotlib.pyplot as plt

from plot_cce_swsh_grid import flatten_grid, make_grid_figure, swsh_grid


def main():
    parser = argparse.ArgumentParser(
        description="Plot the CCE SWSH collocation grid for every field in "
                     "the worldtube HDF5 file, one PNG per field.")
    parser.add_argument("h5file", help="path to the worldtube HDF5 file")
    parser.add_argument("--row", type=int, default=-1,
                         help="timestep row index to plot (default: -1, "
                              "the last row written)")
    parser.add_argument("--l-max", type=int, default=20,
                         help="CCE_LMAX the run used (default 20)")
    parser.add_argument("--radius", type=float, default=None,
                         help="extraction radius, only used to annotate "
                              "titles/axes (plots are always on the unit "
                              "sphere if this is omitted)")
    parser.add_argument("--outdir", default=None,
                         help="output folder for the PNGs (default: "
                              "<h5file-basename>_swsh_plots next to the "
                              "h5 file)")
    args = parser.parse_args()

    theta, phi = swsh_grid(args.l_max)
    theta_flat, phi_flat = flatten_grid(theta, phi)
    n_expected = theta_flat.size

    outdir = args.outdir
    if outdir is None:
        base = os.path.splitext(args.h5file)[0]
        outdir = "%s_swsh_plots" % base
    os.makedirs(outdir, exist_ok=True)

    with h5py.File(args.h5file, "r") as f:
        fields = sorted(name[:-4] for name in f.keys() if name.endswith(".dat"))
        if not fields:
            raise SystemExit("No '*.dat' datasets found in %s" % args.h5file)

        print("Found %d fields; writing plots to %s" % (len(fields), outdir))

        for field in fields:
            ds = f[field + ".dat"]
            row = ds[args.row, :]
            t_val = row[0]
            values = row[1:]

            if values.size != n_expected:
                print("SKIPPING %s: row has %d collocation values but "
                      "--l-max %d implies %d points"
                      % (field, values.size, args.l_max, n_expected))
                continue

            fig = make_grid_figure(theta_flat, phi_flat, values, field,
                                   t_val, args.row, args.l_max, args.radius)
            out_path = os.path.join(outdir, "%s_swsh_grid.png" % field)
            fig.savefig(out_path, dpi=150, bbox_inches="tight")
            plt.close(fig)  # each field opens a new figure -- close it
                             # immediately or memory grows across all 49
            print("  wrote %s" % out_path)

    print("Done: %d PNGs in %s" % (len(fields), outdir))


if __name__ == "__main__":
    main()

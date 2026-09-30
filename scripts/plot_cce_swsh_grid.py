##
# @brief: Plots the SWSH collocation grid from a CCE worldtube HDF5 file
#          (produced by BSSNCtx::writeCceWorldtube()), colored by one field,
#          as a quick visual "proof it's working" check.
#
# Reconstructs the same (theta, phi) grid as
# bssn::cce::generate_swsh_angles()/collocation_offset() in
# BSSN_GR/src/cceWorldtube.cpp: Gauss-Legendre nodes in theta (ascending,
# theta=0 -> pi), equally-spaced phi_j = 2*pi*j/(2*l_max+1), and the
# phi-fastest/theta-slowest flattening offset = phi_index + num_phi*theta_ring.
#
# Usage:
#   python plot_cce_swsh_grid.py <h5file> [--field Lapse] [--row -1]
#                                 [--l-max 20] [--out grid.png]
#
# @date: 2026-09-30

import argparse

import h5py
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)


def num_theta_points(l_max):
    return l_max + 1


def num_phi_points(l_max):
    return 2 * l_max + 1


def swsh_grid(l_max):
    """Returns theta (ascending, 0..pi) and phi (0..2pi, no phase offset),
    matching generate_swsh_angles()'s convention."""
    n_theta = num_theta_points(l_max)
    n_phi = num_phi_points(l_max)
    x, _ = np.polynomial.legendre.leggauss(n_theta)
    theta = np.arccos(np.sort(x)[::-1])  # descending x -> ascending theta
    phi = 2.0 * np.pi * np.arange(n_phi) / n_phi
    return theta, phi


def flatten_grid(theta, phi):
    """Returns (theta_flat, phi_flat), length n_theta*n_phi, in the same
    phi-fastest/theta-slowest order as collocation_offset()."""
    theta_flat = np.repeat(theta, phi.size)
    phi_flat = np.tile(phi, theta.size)
    return theta_flat, phi_flat


def main():
    parser = argparse.ArgumentParser(
        description="Plot the CCE SWSH collocation grid, colored by one field.")
    parser.add_argument("h5file", help="path to the worldtube HDF5 file")
    parser.add_argument("--field", default="Lapse",
                         help="dataset base name to color by (e.g. Lapse, "
                              "gxx, Kxx, Shiftx); default: Lapse")
    parser.add_argument("--row", type=int, default=-1,
                         help="timestep row index to plot (default: -1, "
                              "the last row written)")
    parser.add_argument("--l-max", type=int, default=20,
                         help="CCE_LMAX the run used (default 20)")
    parser.add_argument("--radius", type=float, default=None,
                         help="extraction radius, only used to annotate the "
                              "title (the sphere plot itself is always unit "
                              "radius for clarity)")
    parser.add_argument("--out", default=None,
                         help="output image path (default: "
                              "<field>_swsh_grid.png next to the h5 file)")
    args = parser.parse_args()

    theta, phi = swsh_grid(args.l_max)
    theta_flat, phi_flat = flatten_grid(theta, phi)

    with h5py.File(args.h5file, "r") as f:
        dset_name = args.field + ".dat"
        if dset_name not in f:
            raise SystemExit(
                "Field '%s' not found. Available: %s"
                % (args.field, sorted(n[:-4] for n in f.keys())))
        ds = f[dset_name]
        row = ds[args.row, :]
        t_val = row[0]
        values = row[1:]

    n_expected = theta_flat.size
    if values.size != n_expected:
        raise SystemExit(
            "Row has %d collocation values but --l-max %d implies %d points "
            "-- pass the correct --l-max for this file."
            % (values.size, args.l_max, n_expected))

    # Unit-sphere Cartesian coords for the 3D scatter (radius is cosmetic).
    x = np.sin(theta_flat) * np.cos(phi_flat)
    y = np.sin(theta_flat) * np.sin(phi_flat)
    z = np.cos(theta_flat)

    vmin, vmax = float(values.min()), float(values.max())
    cmap = "viridis"  # perceptually-uniform sequential colormap for magnitude data

    fig = plt.figure(figsize=(13, 6))

    ax3d = fig.add_subplot(1, 2, 1, projection="3d")
    sc = ax3d.scatter(x, y, z, c=values, cmap=cmap, vmin=vmin, vmax=vmax,
                       s=25, depthshade=True)
    ax3d.set_box_aspect((1, 1, 1))
    ax3d.set_xlabel("x")
    ax3d.set_ylabel("y")
    ax3d.set_zlabel("z")
    ax3d.set_title("SWSH collocation grid (unit sphere)")

    ax2d = fig.add_subplot(1, 2, 2)
    sc2 = ax2d.scatter(phi_flat, theta_flat, c=values, cmap=cmap,
                        vmin=vmin, vmax=vmax, s=25)
    ax2d.invert_yaxis()
    ax2d.set_xlabel(r"$\phi$ (equally spaced)")
    ax2d.set_ylabel(r"$\theta$ (Gauss-Legendre, clustered near poles)")
    ax2d.set_title(r"$\theta$-$\phi$ map")

    cbar = fig.colorbar(sc2, ax=[ax3d, ax2d], shrink=0.8, pad=0.05)
    cbar.set_label(args.field)

    radius_note = (" @ r=%g" % args.radius) if args.radius is not None else ""
    fig.suptitle("%s%s, t=%.4g (row %d), l_max=%d (%d points)"
                 % (args.field, radius_note, t_val, args.row, args.l_max,
                    n_expected))

    out = args.out
    if out is None:
        import os
        base = os.path.splitext(args.h5file)[0]
        out = "%s_%s_swsh_grid.png" % (base, args.field)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("Wrote %s" % out)


if __name__ == "__main__":
    main()

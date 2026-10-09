##
# @brief: Injects a known, analytic spherical-harmonic pattern Y_{l,m}(theta,phi)
#          into a copy of an already-validated worldtube HDF5 file, as an
#          empirical check of cceWorldtube.cpp's SWSH angular ordering
#          (Dendro_CCE_v2.1.md Section 9.3 item #2).
#
# Uses m=0 by default so phi-ordering (already confirmed from SpECTRE's own
# SwshCollocation.cpp source -- see project notes) can't interfere: an m=0
# mode has no phi dependence at all, isolating the one still-uncertain
# piece, theta's ring order.
#
# After running this, feed the output through PreprocessCceWorldtube (same
# as any other worldtube file) and check which (l,m) mode actually shows up
# with large amplitude in the resulting modal file (e.g. J.dat's Legend) --
# it should be almost entirely the (l,m) you injected, with other modes at
# the background noise floor. If it shows up at a different l, or spread
# across many modes, that reveals the actual ordering bug.
#
# Usage:
#   python make_ordering_test_worldtube.py <input.h5> <output.h5> [--l 3] [--m 0] [--amplitude 1e-3]
#
# @date: 2026-10-09

import argparse
import shutil

import h5py
import numpy as np
from scipy.special import sph_harm_y

from plot_cce_swsh_grid import swsh_grid, flatten_grid


def main():
    parser = argparse.ArgumentParser(
        description="Inject a known Y_(l,m) pattern into a worldtube file's "
                     "Lapse dataset, for an empirical SWSH-ordering check.")
    parser.add_argument("input_h5", help="path to an existing, validated worldtube h5 file")
    parser.add_argument("output_h5", help="path to write the modified copy")
    parser.add_argument("--l", type=int, default=3,
                         help="l of the injected mode (default 3 -- has 3 "
                              "sign changes in theta, easy to misidentify "
                              "if ring order is wrong)")
    parser.add_argument("--m", type=int, default=0,
                         help="m of the injected mode (default 0, i.e. no "
                              "phi dependence -- isolates theta ordering "
                              "specifically, since phi-ordering is already "
                              "confirmed from SpECTRE's own source)")
    parser.add_argument("--l-max", type=int, default=20,
                         help="CCE_LMAX the input file used (default 20)")
    parser.add_argument("--amplitude", type=float, default=1e-3,
                         help="perturbation amplitude added to Lapse -- "
                              "small relative to Lapse~1, but far above "
                              "the pipeline's existing ~1e-9 noise floor "
                              "so the injected signal clearly dominates")
    args = parser.parse_args()

    theta, phi = swsh_grid(args.l_max)
    theta_flat, phi_flat = flatten_grid(theta, phi)

    # sph_harm_y(l, m, theta, phi): theta=polar/colatitude, phi=azimuthal --
    # SAME convention as this whole project (NOT scipy's older, deprecated
    # sph_harm(m, l, theta, phi), which confusingly swaps theta/phi meaning).
    y_lm = sph_harm_y(args.l, args.m, theta_flat, phi_flat)
    perturbation = args.amplitude * np.real(y_lm)

    shutil.copyfile(args.input_h5, args.output_h5)

    with h5py.File(args.output_h5, "r+") as f:
        lapse = f["Lapse.dat"]
        data = lapse[:, :]
        n_pts = perturbation.size
        if data.shape[1] - 1 != n_pts:
            raise SystemExit(
                "Lapse.dat has %d points but --l-max %d implies %d -- pass "
                "the correct --l-max for this file." % (data.shape[1] - 1,
                                                         args.l_max, n_pts))
        data[:, 1:] += perturbation[np.newaxis, :]
        lapse[...] = data

    print("Wrote %s: injected amplitude=%.3g * Re(Y_%d,%d) into Lapse at "
          "every timestep." % (args.output_h5, args.amplitude, args.l, args.m))
    print("Next: run this through PreprocessCceWorldtube (same as usual), "
          "then check the resulting modal file's per-mode amplitude -- the "
          "(l=%d, m=%d) mode should dominate; if a different l shows up "
          "instead, that's the theta-ring-order bug made visible."
          % (args.l, args.m))


if __name__ == "__main__":
    main()

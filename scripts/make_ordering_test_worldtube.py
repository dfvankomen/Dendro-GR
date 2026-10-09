##
# @brief: Builds a worldtube HDF5 file that is EXACTLY flat space (Lapse=1,
#         g_ii=1, everything else=0) plus one known, analytic
#         spherical-harmonic perturbation Y_{l,m}(theta,phi) added to K's
#         trace (Kxx=Kyy=Kzz), as an empirical check of cceWorldtube.cpp's
#         SWSH angular ordering (Dendro_CCE_v2.1.md Section 9.3 item #2).
#
# Uses m=0 by default so phi-ordering (already confirmed from SpECTRE's own
# SwshCollocation.cpp source -- see project notes) can't interfere: an m=0
# mode has no phi dependence at all, isolating the one still-uncertain
# piece, theta's ring order.
#
# IMPORTANT, three iterations of this script's mistakes so far, all fixed:
# 1. Perturbing the real q1 data gave an inconclusive result -- the real
#    worldtube already has its own physical (l=2,m=+-2) quadrupole content,
#    and the ADM->Bondi gauge transform PreprocessCceWorldtube applies is
#    nonlinear, so a perturbation on top of real data can get mixed with
#    that pre-existing structure, making it impossible to tell which part
#    of the output came from the perturbation vs. the binary's own physics.
#    Fixed by starting from an EXACTLY flat background instead (below).
# 2. Perturbing Lapse's VALUE while leaving DxLapse/DyLapse/DzLapse at their
#    flat value (0) gave an all-zero result -- an angularly-varying Lapse
#    necessarily has a nonzero gradient too, and leaving the derivative
#    fields inconsistent with the value apparently starves the ADM->Bondi
#    conversion of the information it needs. Fixed by switching to K (no
#    paired derivative dataset in the AdmMetricNodal contract -- Section
#    6.1 of Dendro_CCE_v2.1.md -- so nothing to keep consistent).
# 3. Perturbing ONLY Kxx (leaving Kxy/Kyy/Kzz/etc. at zero) gave power
#    spread across many (l,m) with ONLY even m -- because Kxx is a single
#    CARTESIAN TENSOR COMPONENT, not a scalar. A Cartesian component that
#    varies with theta while its tensor partners stay fixed at zero does
#    NOT correspond to an axisymmetric physical field; projecting it onto
#    spherical basis vectors (which themselves carry phi-dependence, e.g.
#    theta_hat ~ (cos(theta)cos(phi), cos(theta)sin(phi), -sin(theta)))
#    mixes in factors like cos^2(phi) -- which decompose into exactly m=0
#    and m=+-2 content, matching what was observed. Fixed by perturbing K's
#    TRACE instead: Kxx=Kyy=Kzz=f(theta,phi) (K_ij = f * delta_ij). Since
#    the identity tensor is the same in any orthonormal basis, this has NO
#    basis-projection mixing -- whatever (l,m) content f has should come
#    through cleanly. Also physically well-motivated: K's trace is tied to
#    the expansion of outgoing null rays, closely related to DuR.
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

from check_cce_worldtube import EXPECTED_DATASETS, GDIAG
from plot_cce_swsh_grid import swsh_grid, flatten_grid

# Flat-space value for each dataset: 1.0 for the metric diagonal and Lapse,
# 0.0 for everything else (off-diagonal metric, all derivatives, Shift, K,
# AuxiliaryShift). Mirrors check_cce_worldtube.py's own near-flat-space
# sanity expectations.
FLAT_VALUE = {name: (1.0 if (name in GDIAG or name == "Lapse") else 0.0)
             for name in EXPECTED_DATASETS}

try:
    from scipy.special import sph_harm_y

    def real_spherical_harmonic(l, m, theta, phi):
        return np.real(sph_harm_y(l, m, theta, phi))
except ImportError:
    # Older scipy (no sph_harm_y yet): fall back to the deprecated
    # sph_harm(m, n, theta, phi), whose argument NAMES are swapped from the
    # physics convention used everywhere else in this project -- scipy's
    # 'theta' is azimuthal (our phi) and its 'phi' is polar/colatitude (our
    # theta). Passing our (theta, phi) straight through would silently
    # compute the WRONG function.
    from scipy.special import sph_harm

    def real_spherical_harmonic(l, m, theta, phi):
        return np.real(sph_harm(m, l, phi, theta))


def main():
    parser = argparse.ArgumentParser(
        description="Inject a known Y_(l,m) pattern into K's trace of a "
                     "worldtube file (over an exactly-flat background), for "
                     "an empirical SWSH-ordering check.")
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
                         help="perturbation amplitude added to K's trace -- "
                              "far above the pipeline's existing ~1e-9 "
                              "noise floor so the injected signal clearly "
                              "dominates")
    args = parser.parse_args()

    theta, phi = swsh_grid(args.l_max)
    theta_flat, phi_flat = flatten_grid(theta, phi)

    # real_spherical_harmonic(l, m, theta, phi): theta=polar/colatitude,
    # phi=azimuthal -- SAME convention as this whole project, regardless of
    # which scipy API is available underneath (see the import block above).
    perturbation = args.amplitude * real_spherical_harmonic(
        args.l, args.m, theta_flat, phi_flat)

    # Only borrowing the input file's structure (time column, row count,
    # Legend, HDF5 attributes) -- every dataset's VALUES get overwritten
    # below with exact flat space, regardless of what the real run had.
    shutil.copyfile(args.input_h5, args.output_h5)

    # Perturb K's TRACE (Kxx=Kyy=Kzz=f), not one lopsided Cartesian
    # component -- see mistake #3 above for why a single component mixes
    # (l,m) content through the basis-vector projection.
    trace_fields = {"Kxx", "Kyy", "Kzz"}

    with h5py.File(args.output_h5, "r+") as f:
        n_pts = perturbation.size
        for name in EXPECTED_DATASETS:
            ds = f[name + ".dat"]
            data = ds[:, :]
            if data.shape[1] - 1 != n_pts:
                raise SystemExit(
                    "%s.dat has %d points but --l-max %d implies %d -- pass "
                    "the correct --l-max for this file."
                    % (name, data.shape[1] - 1, args.l_max, n_pts))
            data[:, 1:] = FLAT_VALUE[name]
            if name in trace_fields:
                data[:, 1:] += perturbation[np.newaxis, :]
            ds[...] = data

    print("Wrote %s: exact flat space + amplitude=%.3g * Re(Y_%d,%d) "
          "injected into K's trace (Kxx=Kyy=Kzz) at every timestep." %
          (args.output_h5, args.amplitude, args.l, args.m))
    print("Next: run this through PreprocessCceWorldtube (same as usual), "
          "then check the resulting modal file's per-mode amplitude -- the "
          "(l=%d, m=%d) mode should dominate; if a different l shows up "
          "instead, that's the theta-ring-order bug made visible."
          % (args.l, args.m))


if __name__ == "__main__":
    main()

##
# @brief: Sanity-checks a CCE worldtube HDF5 file produced by
#          BSSNCtx::writeCceWorldtube() (see Dendro_CCE_v2.1.md Section 9).
# @date: 2026-09-30
#
# Checks, in order:
#   1. The file opens, and contains exactly the 49 expected datasets
#      (structural check against bssn::cce::cce_dataset_names()).
#   2. Every dataset has the expected shape: (num_rows, 1 + N_points),
#      N_points = (l_max+1)*(2*l_max+1) for the run's CCE_LMAX (default 20
#      -> N_points=861, so 862 columns including the time column).
#   3. Column 0 (time) is present, non-decreasing, and matches across all
#      49 datasets row-for-row (they should all be written from the same
#      m_uiTinfo._m_uiT each call).
#   4. No NaN/Inf anywhere in any dataset.
#   5. A physically-motivated sanity check, NOT just "did it run": at
#      r=CCE_EXTRACTION_RADIUS (100M for the default q1 configs) for this
#      near-equal-mass, non-spinning binary, spacetime should be very
#      close to flat (weak-field falloff ~ M/r ~ 0.01 for the metric,
#      ~ M/r^2 ~ 1e-4 for its derivatives). So at every collocation point,
#      every timestep:
#        Lapse            ~ 1       (within LAPSE_TOL)
#        Shift, AuxShift  ~ 0       (within SHIFT_TOL)
#        g_ii (diagonal)  ~ 1       (within GDIAG_TOL)
#        g_ij (off-diag)  ~ 0       (within GOFFDIAG_TOL)
#        K_ij             ~ 0       (within K_TOL)
#        all derivative datasets (Dx/Dy/Dz of g and Shift, DxLapse etc.) ~ 0
#   This will NOT catch every possible bug (e.g. it can't catch the SWSH
#   angular-ordering question flagged in Dendro_CCE_v2.1.md Section 9.3
#   item 2, since flat space has no angular structure to get the ordering
#   wrong about) -- it's a first "is the pipeline basically sane" gate,
#   not a substitute for the analytic test suite in Section 8.

import argparse
import math
import sys

import h5py
import numpy as np

# mirrors bssn::cce::cce_dataset_names() in BSSN_GR/src/cceWorldtube.cpp
EXPECTED_DATASETS = [
    "gxx", "gxy", "gxz", "gyy", "gyz", "gzz",
    "Dxgxx", "Dxgxy", "Dxgxz", "Dxgyy", "Dxgyz", "Dxgzz",
    "Dygxx", "Dygxy", "Dygxz", "Dygyy", "Dygyz", "Dygzz",
    "Dzgxx", "Dzgxy", "Dzgxz", "Dzgyy", "Dzgyz", "Dzgzz",
    "Shiftx", "Shifty", "Shiftz",
    "DxShiftx", "DxShifty", "DxShiftz",
    "DyShiftx", "DyShifty", "DyShiftz",
    "DzShiftx", "DzShifty", "DzShiftz",
    "Lapse",
    "DxLapse", "DyLapse", "DzLapse",
    "Kxx", "Kxy", "Kxz", "Kyy", "Kyz", "Kzz",
    "AuxiliaryShiftx", "AuxiliaryShifty", "AuxiliaryShiftz",
]

GDIAG = {"gxx", "gyy", "gzz"}
GOFFDIAG = {"gxy", "gxz", "gyz"}
SHIFT_LIKE = {"Shiftx", "Shifty", "Shiftz",
             "AuxiliaryShiftx", "AuxiliaryShifty", "AuxiliaryShiftz"}
K_LIKE = {"Kxx", "Kxy", "Kxz", "Kyy", "Kyz", "Kzz"}
DERIV_LIKE = {name for name in EXPECTED_DATASETS if name.startswith(("Dx", "Dy", "Dz"))}


def num_collocation_points(l_max):
    return (l_max + 1) * (2 * l_max + 1)


def main():
    parser = argparse.ArgumentParser(
        description="Sanity-check a CCE worldtube HDF5 file.")
    parser.add_argument("h5file", help="path to the worldtube HDF5 file "
                                        "(e.g. cce_worldtube_test.h5)")
    parser.add_argument("--l-max", type=int, default=20,
                         help="CCE_LMAX the run used (default 20, matching "
                              "the compiled-in default).")
    parser.add_argument("--lapse-tol", type=float, default=0.05)
    parser.add_argument("--shift-tol", type=float, default=0.01)
    parser.add_argument("--gdiag-tol", type=float, default=0.05)
    parser.add_argument("--goffdiag-tol", type=float, default=0.02)
    parser.add_argument("--k-tol", type=float, default=0.01)
    parser.add_argument("--deriv-tol", type=float, default=0.01)
    args = parser.parse_args()

    expected_n_pts = num_collocation_points(args.l_max)
    expected_ncols = expected_n_pts + 1

    failures = []
    warnings = []

    print("Opening %s ..." % args.h5file)
    f = h5py.File(args.h5file, "r")

    # --- 1. structural: dataset presence ---
    present = set(f.keys())
    expected_names = {name + ".dat" for name in EXPECTED_DATASETS}
    missing = expected_names - present
    extra = present - expected_names
    if missing:
        failures.append("Missing datasets: %s" % sorted(missing))
    if extra:
        warnings.append("Unexpected extra datasets present: %s" % sorted(extra))
    print("Dataset presence: %d/%d expected datasets found%s"
          % (len(expected_names) - len(missing), len(expected_names),
             "" if not missing else " (MISSING SOME)"))

    if missing:
        print("\nFAIL: cannot continue, missing datasets.")
        for msg in failures:
            print("  - " + msg)
        sys.exit(1)

    # --- 2. shape check ---
    n_rows = None
    for name in EXPECTED_DATASETS:
        ds = f[name + ".dat"]
        shape = ds.shape
        if shape[1] != expected_ncols:
            failures.append(
                "%s.dat has %d columns, expected %d (= 1 time + %d "
                "collocation points at l_max=%d)"
                % (name, shape[1], expected_ncols, expected_n_pts, args.l_max))
        if n_rows is None:
            n_rows = shape[0]
        elif shape[0] != n_rows:
            failures.append(
                "%s.dat has %d rows, but gxx.dat (checked first) has %d -- "
                "row counts should match across all 49 datasets"
                % (name, shape[0], n_rows))
    print("Row count: %d timesteps written (per-dataset)" % (n_rows or 0))
    if not n_rows:
        print("\nFAIL: zero rows written -- CCE_ENABLED may not have taken "
              "effect, or the run didn't reach CCE_OUTPUT_FREQ steps.")
        sys.exit(1)

    # --- 3. time column consistency ---
    time_ref = f["gxx.dat"][:, 0]
    if np.any(np.diff(time_ref) < 0):
        failures.append("Time column in gxx.dat is not non-decreasing")
    for name in EXPECTED_DATASETS[1:]:
        t = f[name + ".dat"][:, 0]
        if not np.allclose(t, time_ref):
            failures.append(
                "%s.dat's time column doesn't match gxx.dat's row-for-row"
                % name)
    print("Time range: [%g, %g] over %d rows"
          % (time_ref[0], time_ref[-1], len(time_ref)))

    # --- 4 & 5: NaN/Inf and physical sanity, per dataset ---
    print("\nPer-dataset checks (expected near-flat-space values at "
          "r=%s extraction radius):" % "CCE_EXTRACTION_RADIUS")
    print("%-20s %10s %14s %14s %10s" %
          ("dataset", "expected~", "min", "max", "status"))

    for name in EXPECTED_DATASETS:
        data = f[name + ".dat"][:, 1:]  # drop time column
        finite = np.isfinite(data)
        if not finite.all():
            n_bad = (~finite).sum()
            failures.append("%s.dat has %d non-finite (NaN/Inf) values"
                             % (name, n_bad))
            print("%-20s %10s %14s %14s %10s" %
                  (name, "-", "-", "-", "NaN/Inf!"))
            continue

        if name in GDIAG:
            expected, tol = 1.0, args.gdiag_tol
        elif name in GOFFDIAG:
            expected, tol = 0.0, args.goffdiag_tol
        elif name == "Lapse":
            expected, tol = 1.0, args.lapse_tol
        elif name in SHIFT_LIKE:
            expected, tol = 0.0, args.shift_tol
        elif name in K_LIKE:
            expected, tol = 0.0, args.k_tol
        elif name in DERIV_LIKE:
            expected, tol = 0.0, args.deriv_tol
        else:
            expected, tol = None, None

        dmin, dmax = float(data.min()), float(data.max())
        status = "ok"
        if expected is not None:
            if dmin < expected - tol or dmax > expected + tol:
                status = "OUT OF RANGE"
                failures.append(
                    "%s.dat: values in [%.6g, %.6g], expected ~%.3g +/- %.3g "
                    "(near-flat-space sanity check)"
                    % (name, dmin, dmax, expected, tol))
        print("%-20s %10s %14.6g %14.6g %10s" %
              (name, ("%.2g" % expected) if expected is not None else "-",
               dmin, dmax, status))

    f.close()

    print()
    if failures:
        print("FAIL: %d issue(s) found:" % len(failures))
        for msg in failures:
            print("  - " + msg)
        if warnings:
            print("Also, %d warning(s):" % len(warnings))
            for msg in warnings:
                print("  - " + msg)
        sys.exit(1)
    else:
        print("PASS: structure, time consistency, finiteness, and "
              "near-flat-space sanity checks all passed.")
        if warnings:
            print("(%d non-fatal warning(s):" % len(warnings))
            for msg in warnings:
                print("  - " + msg)
        print("\nNote: this does NOT confirm the file is byte-correct for "
              "SpECTRE's PreprocessCceWorldtube to read (dataset attribute "
              "format, HDF5 Legend/version attributes, and the SWSH angular "
              "ordering are not checked here) -- see Dendro_CCE_v2.1.md "
              "Section 9.3 items 2 and 7 for what's still unverified.")
        sys.exit(0)


if __name__ == "__main__":
    main()

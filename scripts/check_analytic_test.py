##
# @brief: Compares AnalyticTestCharacteristicExtract's actual News output
#          against the exact analytic answer bundled in the same file
#          (Cce/News_expected.dat), using the tolerance declared in the
#          test's own YAML (AnalyticTestLinearizedBondiSachs.yaml's
#          OutputFileChecks: AbsoluteTolerance: 1e-2).
#
# Usage:
#   python check_analytic_test.py CharacteristicExtractReduction.h5
#
# @date: 2026-10-09

import argparse

import h5py
import numpy as np


def main():
    parser = argparse.ArgumentParser(
        description="Compare an AnalyticTestCharacteristicExtract run's "
                     "News against its bundled exact analytic answer.")
    parser.add_argument("h5file", help="path to CharacteristicExtractReduction.h5")
    parser.add_argument("--abs-tol", type=float, default=1e-2,
                         help="absolute tolerance (default 1e-2, matching "
                              "AnalyticTestLinearizedBondiSachs.yaml's own "
                              "OutputFileChecks declaration)")
    args = parser.parse_args()

    with h5py.File(args.h5file, "r") as f:
        cce_group = [k for k in f.keys() if k.endswith(".cce")][0]
        actual = f[f"{cce_group}/News"][:, :]
        expected = f["Cce/News_expected.dat"][:, :]

    if actual.shape != expected.shape:
        raise SystemExit("Shape mismatch: actual %s vs expected %s -- "
                         "can't compare directly." % (actual.shape, expected.shape))

    t_actual, t_expected = actual[:, 0], expected[:, 0]
    if not np.allclose(t_actual, t_expected):
        raise SystemExit("Time columns don't match between actual and "
                         "expected -- rows aren't aligned, comparison would "
                         "be meaningless.")

    diff = np.abs(actual[:, 1:] - expected[:, 1:])
    max_diff = diff.max()
    worst_row, worst_col = np.unravel_index(diff.argmax(), diff.shape)

    print("Rows compared: %d, time range [%.4g, %.4g]"
          % (actual.shape[0], t_actual.min(), t_actual.max()))
    print("Max |actual - expected| = %.6g at row %d (t=%.4g), mode column %d"
          % (max_diff, worst_row, t_actual[worst_row], worst_col + 1))
    print("Tolerance: %.6g" % args.abs_tol)

    if max_diff <= args.abs_tol:
        print("PASS: within tolerance.")
    else:
        print("FAIL: exceeds tolerance.")


if __name__ == "__main__":
    main()

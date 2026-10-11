##
# @brief: Repairs a worldtube HDF5 file that was copied while Dendro-GR was
#          still actively writing to it (a torn-write race, not a code bug --
#          writeCceWorldtube() extends each dataset's size before writing its
#          data and only flushes once per timestep after all 49 datasets are
#          appended, so a copy landing mid-timestep can capture a row whose
#          size was extended but whose data was never actually written,
#          reading back as zero -- including a zeroed time column, which
#          SpECTRE's WorldtubeBufferUpdater correctly rejects as
#          non-monotonic).
#
# Finds, for every one of the 49 expected datasets, the first row index
# where its time column stops being strictly increasing (if any), takes the
# minimum such safe row count across all of them, and writes a new file
# truncated uniformly to that many rows -- so every dataset stays
# row-for-row time-aligned, matching what check_cce_worldtube.py itself
# verifies.
#
# Usage:
#   python fix_torn_write.py <input.h5> <output.h5>
#
# @date: 2026-10-10

import argparse
import shutil

import h5py
import numpy as np

from check_cce_worldtube import EXPECTED_DATASETS


def first_non_monotonic_index(t):
    """Returns the row index of the first break in strict monotonicity, or
    len(t) if the whole column is clean."""
    diffs = np.diff(t)
    bad = np.where(diffs <= 0)[0]
    return int(bad[0] + 1) if bad.size else len(t)


def main():
    parser = argparse.ArgumentParser(
        description="Truncate a torn-write-corrupted worldtube h5 file to "
                     "its last fully-consistent row, uniformly across all "
                     "49 datasets.")
    parser.add_argument("input_h5")
    parser.add_argument("output_h5")
    args = parser.parse_args()

    safe_rows = None
    with h5py.File(args.input_h5, "r") as f:
        for name in EXPECTED_DATASETS:
            t = f[name + ".dat"][:, 0]
            n = first_non_monotonic_index(t)
            if n < len(t):
                print("%s.dat: time column breaks monotonicity at row %d "
                      "(t=%.6g -> %.6g)" % (name, n, t[n - 1], t[n]))
            if safe_rows is None or n < safe_rows:
                safe_rows = n

    print("Safe row count across all datasets: %d" % safe_rows)

    shutil.copyfile(args.input_h5, args.output_h5)
    with h5py.File(args.output_h5, "r+") as f:
        for name in EXPECTED_DATASETS:
            ds = f[name + ".dat"]
            data = ds[:safe_rows, :]
            ds.resize((safe_rows, ds.shape[1]))
            ds[...] = data

    print("Wrote %s, truncated to %d rows (t up to %.6g)"
         % (args.output_h5, safe_rows, data[-1, 0]))


if __name__ == "__main__":
    main()

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
# Unlike the first version of this script (which just truncated everything
# after the FIRST bad row), this scans the WHOLE file and removes every bad
# row it finds, wherever it occurs -- in case a long multi-day run gets
# copied mid-write more than once. For each of the 49 expected datasets, a
# row is "bad" if its time value doesn't exceed the time of the last row
# RETAINED so far (checked per-dataset, against the retained sequence, so
# multiple scattered incidents are each caught correctly rather than only
# the first). The union of bad row indices across all 49 datasets is then
# dropped uniformly from every dataset, so they all stay row-for-row
# time-aligned, matching what check_cce_worldtube.py itself verifies.
#
# Usage:
#   python fix_torn_write.py <input.h5> <output.h5>
#
# @date: 2026-10-13

import argparse
import shutil

import h5py
import numpy as np

from check_cce_worldtube import EXPECTED_DATASETS


def find_bad_rows(t):
    """Returns the set of row indices in time column `t` that are bad --
    non-increasing relative to the last RETAINED good row, not just the
    immediately preceding row, so consecutive/scattered bad rows are each
    caught rather than only the first one found."""
    bad = set()
    last_good_time = t[0]
    for i in range(1, len(t)):
        if t[i] > last_good_time:
            last_good_time = t[i]
        else:
            bad.add(i)
    return bad


def main():
    parser = argparse.ArgumentParser(
        description="Remove every torn-write-corrupted row from a worldtube "
                     "h5 file, wherever it occurs, uniformly across all 49 "
                     "datasets.")
    parser.add_argument("input_h5")
    parser.add_argument("output_h5")
    args = parser.parse_args()

    bad_rows = set()
    n_rows = None
    with h5py.File(args.input_h5, "r") as f:
        for name in EXPECTED_DATASETS:
            t = f[name + ".dat"][:, 0]
            n_rows = len(t) if n_rows is None else n_rows
            this_bad = find_bad_rows(t)
            for i in sorted(this_bad):
                print("%s.dat: row %d is bad (t=%.6g, following t=%.6g)"
                     % (name, i, t[i], t[i - 1]))
            bad_rows |= this_bad

    if not bad_rows:
        print("No bad rows found -- file is already clean.")
    else:
        print("Removing %d bad row(s) total (union across all datasets): %s"
             % (len(bad_rows), sorted(bad_rows)))

    keep = np.array([i for i in range(n_rows) if i not in bad_rows])

    shutil.copyfile(args.input_h5, args.output_h5)
    with h5py.File(args.output_h5, "r+") as f:
        for name in EXPECTED_DATASETS:
            ds = f[name + ".dat"]
            data = ds[:, :][keep, :]
            ds.resize((len(keep), ds.shape[1]))
            ds[...] = data

    print("Wrote %s: %d rows kept (of %d), t up to %.6g"
         % (args.output_h5, len(keep), n_rows, data[-1, 0]))


if __name__ == "__main__":
    main()

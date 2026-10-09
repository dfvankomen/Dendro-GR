##
# @brief: Plots the same (l,m) mode for EVERY field in a CharacteristicExtract
#          output file's '*.cce' group (News, Strain, Psi0-4,
#          EthInertialRetardedTime, ...), saving one PNG per field into an
#          output folder. Same plot (Re/Im/norm vs retarded time) as
#          plot_cce_news.py, just looped over every field instead of picking
#          one with --field.
#
# Usage:
#   python plotall_cce_fields.py <h5file> [--l 2] [--m 2] [--outdir cce_plots]
#
# @date: 2026-10-09

import argparse
import os

import h5py
import matplotlib
matplotlib.use("Agg")  # no display needed; also avoids piling up open GUI
                        # windows across every field
import matplotlib.pyplot as plt

from plot_cce_news import find_cce_group, legend_column, make_mode_figure


def main():
    parser = argparse.ArgumentParser(
        description="Plot the same (l,m) mode for every field in a "
                     "CharacteristicExtract output file, one PNG per field.")
    parser.add_argument("h5file", help="path to CharacteristicExtractReduction.h5")
    parser.add_argument("--l", type=int, default=2, help="l mode (default 2)")
    parser.add_argument("--m", type=int, default=2, help="m mode (default 2)")
    parser.add_argument("--outdir", default=None,
                         help="output folder for the PNGs (default: "
                              "<h5file-basename>_cce_plots next to the "
                              "h5 file)")
    args = parser.parse_args()

    outdir = args.outdir
    if outdir is None:
        base = os.path.splitext(args.h5file)[0]
        outdir = "%s_cce_plots" % base
    os.makedirs(outdir, exist_ok=True)

    with h5py.File(args.h5file, "r") as f:
        cce_group = find_cce_group(f)
        fields = sorted(name for name, obj in f[cce_group].items()
                        if isinstance(obj, h5py.Dataset))

    print("Found %d fields in '%s'; writing plots to %s"
          % (len(fields), cce_group, outdir))

    written = 0
    for field in fields:
        with h5py.File(args.h5file, "r") as f:
            ds = f["%s/%s" % (cce_group, field)]
            legend = ds.attrs.get("Legend")
            if legend is None:
                print("SKIPPING %s: no Legend attribute" % field)
                continue
            try:
                re_col = legend_column(legend, "Real", args.l, args.m)
                im_col = legend_column(legend, "Imag", args.l, args.m)
            except SystemExit:
                print("SKIPPING %s: no (l=%d, m=%d) mode found"
                      % (field, args.l, args.m))
                continue
            data = ds[:, :]

        t = data[:, 0]
        re = data[:, re_col]
        im = data[:, im_col]

        fig = make_mode_figure(t, re, im, field, args.l, args.m)
        out_path = os.path.join(outdir, "%s_l%dm%d.png" % (field, args.l, args.m))
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)  # each field opens a new figure -- close it
                         # immediately or memory grows across all fields
        print("  wrote %s" % out_path)
        written += 1

    print("Done: %d/%d PNGs in %s" % (written, len(fields), outdir))


if __name__ == "__main__":
    main()

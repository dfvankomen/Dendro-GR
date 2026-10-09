##
# @brief: Plots one (l,m) mode of the Bondi News (or Strain/Psi0-4) from a
#          CharacteristicExtract output file, as a "proof the CCE pipeline
#          produces a real waveform" check. No BMS/scri frame-fixing -- this
#          is the raw CCE output, good enough to show the pipeline works,
#          not yet science-ready for comparison across frames/gauges.
#
# Usage:
#   python plot_cce_news.py <h5file> [--field News] [--l 2] [--m 2] [--out news_l2m2.png]
#
# @date: 2026-10-09

import argparse
import re

import h5py
import matplotlib.pyplot as plt
import numpy as np


def find_cce_group(f):
    """Finds the single '<Name>.cce' group in the file (holds News, Strain,
    Psi0-4, EthInertialRetardedTime) rather than hardcoding a radius-specific
    name like 'SpectreR0100.cce'."""
    candidates = [name for name in f.keys() if name.endswith(".cce")]
    if len(candidates) != 1:
        raise SystemExit(
            "Expected exactly one '*.cce' group, found: %s" % candidates)
    return candidates[0]


def legend_column(legend, prefix, l, m):
    """Finds the column index in `legend` matching '<prefix> Y_l,m', e.g.
    legend_column(legend, 'Real', 2, 2) -> index of 'Real Y_2,2'. This is
    CharacteristicExtract's own convention (confirmed from a real Legend
    attribute) -- NOT the same 'Re(l,m)'/'Im(l,m)' convention used by
    PreprocessCceWorldtube's output earlier in the pipeline."""
    target = "%s Y_%d,%d" % (prefix, l, m)
    for i, name in enumerate(legend):
        name_str = name.decode() if isinstance(name, bytes) else name
        if name_str == target:
            return i
    raise SystemExit(
        "Column '%s' not found in Legend. Available: %s"
        % (target, [n.decode() if isinstance(n, bytes) else n
                    for n in legend]))


def read_mode(h5file, field, l, m):
    """Returns (t, re, im) for one (l,m) mode of `field` in the single
    '*.cce' group of `h5file`. Shared by plot_cce_news.py (one field) and
    plotall_cce_fields.py (loops over every field) so they can't drift out
    of sync."""
    with h5py.File(h5file, "r") as f:
        cce_group = find_cce_group(f)
        ds = f["%s/%s" % (cce_group, field)]
        legend = ds.attrs["Legend"]
        data = ds[:, :]

    re_col = legend_column(legend, "Real", l, m)
    im_col = legend_column(legend, "Imag", l, m)
    return data[:, 0], data[:, re_col], data[:, im_col]


def make_mode_figure(t, re, im, field, l, m):
    """Builds and returns the Re/Im/norm-vs-time figure for one (l,m) mode
    of one CCE output quantity."""
    norm = np.sqrt(re**2 + im**2)

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(t, re, color="#4C72B0", marker="o", markersize=4,
           label="Re")
    ax.plot(t, im, color="#DD8452", marker="o", markersize=4,
           label="Im")
    ax.plot(t, norm, color="#55A868", linestyle="--", marker="o",
           markersize=4, label="norm")
    ax.set_xlabel("retarded time (M)")
    ax.set_ylabel("%s, (l=%d, m=%d) mode" % (field, l, m))
    ax.set_title("%s (l=%d, m=%d) at scri+ -- raw CCE output, no BMS frame fix"
                % (field, l, m))
    ax.legend()
    ax.grid(alpha=0.3)
    return fig


def main():
    parser = argparse.ArgumentParser(
        description="Plot one (l,m) mode of a CCE waveform quantity "
                     "(News, Strain, Psi0-4) vs time.")
    parser.add_argument("h5file", help="path to CharacteristicExtractReduction.h5")
    parser.add_argument("--field", default="News",
                         help="which quantity to plot (News, Strain, Psi0, "
                              "Psi1, Psi2, Psi3, Psi4); default: News")
    parser.add_argument("--l", type=int, default=2, help="l mode (default 2)")
    parser.add_argument("--m", type=int, default=2, help="m mode (default 2)")
    parser.add_argument("--out", default=None,
                         help="output image path (default: "
                              "<field>_l<l>m<m>.png next to the h5 file)")
    args = parser.parse_args()

    t, re, im = read_mode(args.h5file, args.field, args.l, args.m)
    fig = make_mode_figure(t, re, im, args.field, args.l, args.m)

    out = args.out
    if out is None:
        import os
        base = os.path.splitext(args.h5file)[0]
        out = "%s_%s_l%dm%d.png" % (base, args.field, args.l, args.m)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("Wrote %s" % out)


if __name__ == "__main__":
    main()

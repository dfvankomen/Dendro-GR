"""Compare constraint norms of outer-BC runs (dgr_ConstraintNorms.dat, plain or .gz).

usage: python3 compare.py label=path [label=path ...] [--col C1,C2] [--png out.png]

Prints each column at t = 1, 5, 10, 20, 30, ... M for every run, then the ratio of every
run to the first. Default columns: C_HAM_rms_far, C_MOM_MAG_rms_far (r > BSSN_CONSTRAINT_FAR_R).
  python3 compare.py mask5tail=ref/A6_mask5tail.dat.gz koskip4=ref/A6_koskip4.dat.gz \\
      ghost=<run dir>/ghost_A6/dat/dgr_ConstraintNorms.dat --png A6.png
"""
import argparse
import gzip

import numpy as np

REPORT_TIMES = [1, 5, 10, 20, 30, 40, 50, 60, 80, 100, 150, 200, 300, 400]


def load_norms(path):
    """Columns of a ConstraintNorms file as a dict of arrays, keyed by header name."""
    with (gzip.open if path.endswith(".gz") else open)(path, "rt") as f:
        rows = [line.split() for line in f if line.strip() and not line.startswith("#")]
    data = np.array(rows[1:], dtype=float)
    return {name: data[:, i] for i, name in enumerate(rows[0])}


def value_at(run, column, t):
    """Column interpolated to time t; NaN if the run ended more than 0.5 M before t."""
    t_end = run["time"].max()
    return np.interp(min(t, t_end), run["time"], run[column]) if t_end >= t - 0.5 else np.nan


def print_table(runs, column):
    labels = list(runs)
    base = labels[0]
    print(f"\n{column}")
    print("  t/M" + "".join(f"{l:>14}" for l in labels) + "".join(f"  {l + '/' + base:>22}" for l in labels[1:]))
    t_last = max(r["time"].max() for r in runs.values())
    for t in (t for t in REPORT_TIMES if t <= t_last + 0.5):
        values = [value_at(runs[l], column, t) for l in labels]
        ratios = [v / values[0] for v in values[1:]]
        print(f"{t:5.0f}" + "".join(f"{v:14.3e}" for v in values) + "".join(f"  {r:22.3g}" for r in ratios))


def plot(runs, columns, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, len(columns), figsize=(5.5 * len(columns), 4), squeeze=False)
    for ax, column in zip(axes[0], columns):
        for label, run in runs.items():
            ax.semilogy(run["time"], run[column], lw=1.5, label=label)
        ax.set_xlabel("t / M")
        ax.set_title(column, fontsize=10)
        ax.grid(alpha=0.3, which="both", lw=0.5)
    axes[0][0].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    print(f"\nwrote {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("runs", nargs="+", metavar="label=path")
    parser.add_argument("--col", default="C_HAM_rms_far,C_MOM_MAG_rms_far")
    parser.add_argument("--png")
    args = parser.parse_args()

    runs = {}
    for spec in args.runs:
        label, path = spec.split("=", 1)
        runs[label] = load_norms(path)
    columns = args.col.split(",")
    for column in columns:
        print_table(runs, column)
    if args.png:
        plot(runs, columns, args.png)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Pick the VTU files from a Dendro-GR run that fall inside a time window.

Each time Dendro-GR writes VTU output it makes one .pvtu file (the one you open
in ParaView) and one .vtu file per MPI rank:

    bssn_gr_122500.pvtu
    bssn_gr_122500_0_256.vtu ... bssn_gr_122500_255_256.vtu

The number after the prefix is the time step. This script reads dgr_GridInfo.dat
to turn time steps into simulation times, keeps the steps between --tmin and
--tmax, and then shows, lists, copies, links or tars their files.

Examples, run from a run directory that has vtu/ and dat/ in it:

  # see which steps would be picked, with times and sizes
  select_vtu.py --vtu-dir vtu --tmin 600 --tmax 850

  # write the file list, then copy just those files to your own machine
  select_vtu.py --vtu-dir vtu --tmin 600 --tmax 850 --mode list > files.txt
  rsync -av --files-from=user@cluster:/path/run/files.txt user@cluster:/path/run/vtu/ ./late_vtu/

  # put every 2nd step into one tar file
  select_vtu.py --vtu-dir vtu --tmin 600 --tmax 850 --every 2 --mode tar --dest late.tar
"""
import argparse
import bisect
import glob
import os
import re
import shutil
import sys
import tarfile


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vtu-dir", required=True,
                   help="folder with the .pvtu and .vtu files")
    p.add_argument("--dat", default=None,
                   help="file that maps time steps to times (default: dgr_GridInfo.dat, "
                        "looked for in VTU_DIR, VTU_DIR/../dat, ./dat and .)")
    p.add_argument("--tmin", type=float, default=float("-inf"),
                   help="earliest time to keep, in M (default: start of the run)")
    p.add_argument("--tmax", type=float, default=float("inf"),
                   help="latest time to keep, in M (default: end of the run)")
    p.add_argument("--every", type=int, default=1,
                   help="keep only every Nth picked step, e.g. 2 halves the data (default: 1)")
    p.add_argument("--mode", choices=["show", "list", "copy", "link", "tar"], default="show",
                   help="what to do with the picked files. "
                        "show: print a table of steps, times and sizes. "
                        "list: print file names, one per line (for rsync --files-from). "
                        "copy / link: copy or symlink them into the folder DEST. "
                        "tar: pack them into the tar file DEST. (default: show)")
    p.add_argument("--dest", default=None,
                   help="output folder for copy/link, or output file for tar")
    p.add_argument("--no-pvtu", action="store_true",
                   help="leave out the .pvtu files and keep only the per-rank .vtu files")

    names = p.add_argument_group("file naming (only needed if your run names things differently)")
    names.add_argument("--prefix", default=None,
                       help="file name prefix, i.e. the end of BSSN_VTU_FILE_PREFIX "
                            "(default: worked out from the .pvtu files)")
    names.add_argument("--pvtu-pattern", default="{prefix}_{step}.pvtu",
                       help="name of the .pvtu files (default: %(default)s)")
    names.add_argument("--vtu-pattern", default="{prefix}_{step}_*.vtu",
                       help="name of the per-rank .vtu files, * matches anything (default: %(default)s)")
    names.add_argument("--step-col", default="timeStep",
                       help="name of the time-step column in DAT (default: %(default)s)")
    names.add_argument("--time-col", default="simTime",
                       help="name of the time column in DAT (default: %(default)s)")
    return p.parse_args()


def find_dat(vtu_dir):
    """Look for dgr_GridInfo.dat in the usual places."""
    for folder in (vtu_dir, os.path.join(vtu_dir, "..", "dat"), "dat", "."):
        path = os.path.join(folder, "dgr_GridInfo.dat")
        if os.path.isfile(path):
            return os.path.normpath(path)
    sys.exit(f"error: could not find dgr_GridInfo.dat near {vtu_dir}; pass it with --dat")


def read_times(path, step_col, time_col):
    """Read (step, time) pairs from a CSV file; return sorted steps and their times.

    After a restart from a checkpoint the file repeats its header line and some
    steps. A repeated step keeps the time from its last row.
    """
    time_at_step = {}
    columns = None
    with open(path) as f:
        for line in f:
            fields = [x.strip() for x in line.split(",")]
            if step_col in fields:  # a header line
                if time_col not in fields:
                    sys.exit(f"error: {path} has no '{time_col}' column")
                columns = (fields.index(step_col), fields.index(time_col))
                continue
            if columns is None:
                continue
            try:
                time_at_step[int(float(fields[columns[0]]))] = float(fields[columns[1]])
            except (ValueError, IndexError):
                continue  # blank or broken line
    if not time_at_step:
        sys.exit(f"error: no '{step_col}' / '{time_col}' data found in {path}")
    steps = sorted(time_at_step)
    return steps, [time_at_step[s] for s in steps]


def time_of(step, steps, times):
    """Time of a step: (time, True) if the file lists it, (estimate, False) if it lies
    between two listed steps, (None, False) if it is outside the file's range."""
    i = bisect.bisect_left(steps, step)
    if i < len(steps) and steps[i] == step:
        return times[i], True
    if i == 0 or i == len(steps):
        return None, False
    s0, s1, t0, t1 = steps[i - 1], steps[i], times[i - 1], times[i]
    return t0 + (t1 - t0) * (step - s0) / (s1 - s0), False


def find_output_steps(vtu_dir, pvtu_pattern, prefix):
    """Return the file prefix and the sorted time steps that have a .pvtu file."""
    # turn e.g. "{prefix}_{step}.pvtu" into a regular expression
    regex = re.escape(pvtu_pattern)
    regex = regex.replace(r"\{prefix\}", r"(?P<prefix>.+?)").replace(r"\{step\}", r"(?P<step>\d+)")
    regex = re.compile("^" + regex + "$")

    steps_by_prefix = {}
    for name in os.listdir(vtu_dir):
        m = regex.match(name)
        if m and (prefix is None or m.group("prefix") == prefix):
            steps_by_prefix.setdefault(m.group("prefix"), []).append(int(m.group("step")))

    if not steps_by_prefix:
        sys.exit(f"error: no files named like '{pvtu_pattern}' in {vtu_dir}")
    if len(steps_by_prefix) > 1:
        found = ", ".join(sorted(steps_by_prefix))
        sys.exit(f"error: found more than one prefix ({found}); choose one with --prefix")
    prefix, steps = next(iter(steps_by_prefix.items()))
    return prefix, sorted(steps)


def main():
    args = parse_args()
    if args.every < 1:
        sys.exit("error: --every must be 1 or more")
    if args.mode in ("copy", "link", "tar") and not args.dest:
        sys.exit(f"error: --mode {args.mode} needs --dest")
    dat = args.dat or find_dat(args.vtu_dir)

    prefix, output_steps = find_output_steps(args.vtu_dir, args.pvtu_pattern, args.prefix)
    steps, times = read_times(dat, args.step_col, args.time_col)

    # keep the output steps whose time is inside [tmin, tmax]
    picked = []
    n_outside = 0
    for step in output_steps:
        t, exact = time_of(step, steps, times)
        if t is None:
            n_outside += 1
        elif args.tmin <= t <= args.tmax:
            picked.append((step, t, exact))
    picked = picked[::args.every]
    if not picked:
        sys.exit(f"error: no output steps between t = {args.tmin} and t = {args.tmax}")
    n_estimated = sum(not exact for _, _, exact in picked)

    # collect the files of each picked step
    selection = []  # (step, time, files, bytes)
    for step, t, _ in picked:
        pattern = args.vtu_pattern.format(prefix=prefix, step=step)
        files = sorted(glob.glob(os.path.join(args.vtu_dir, pattern)))
        if not args.no_pvtu:
            files.insert(0, os.path.join(args.vtu_dir,
                                         args.pvtu_pattern.format(prefix=prefix, step=step)))
        selection.append((step, t, files, sum(os.path.getsize(f) for f in files)))
    n_files = sum(len(files) for _, _, files, _ in selection)
    gigabytes = sum(size for _, _, _, size in selection) / 1e9

    if args.mode == "show":
        print(f"{prefix}: {len(selection)} of {len(output_steps)} output steps picked, "
              f"t = {selection[0][1]:.2f} to {selection[-1][1]:.2f}, "
              f"{n_files} files, {gigabytes:.2f} GB (times from {dat})")
        print(f"{'step':>8} {'t [M]':>10} {'files':>6} {'MB':>9}")
        for step, t, files, size in selection:
            print(f"{step:8d} {t:10.3f} {len(files):6d} {size / 1e6:9.1f}")

    elif args.mode == "list":
        for _, _, files, _ in selection:
            for f in files:
                print(os.path.relpath(f, args.vtu_dir))

    elif args.mode in ("copy", "link"):
        os.makedirs(args.dest, exist_ok=True)
        for _, _, files, _ in selection:
            for f in files:
                target = os.path.join(args.dest, os.path.basename(f))
                if args.mode == "copy":
                    shutil.copy2(f, target)
                elif not os.path.lexists(target):
                    os.symlink(os.path.abspath(f), target)
        print(f"{args.mode}: {n_files} files ({gigabytes:.2f} GB) -> {args.dest}", file=sys.stderr)

    else:  # tar
        with tarfile.open(args.dest, "w") as tar:
            for _, _, files, _ in selection:
                for f in files:
                    tar.add(f, arcname=os.path.basename(f))
        print(f"tar: {n_files} files ({gigabytes:.2f} GB) -> {args.dest}", file=sys.stderr)

    # notes go to stderr so they don't end up in a --mode list file
    if n_estimated:
        print(f"note: {n_estimated} picked step(s) are not listed in {dat}; "
              "their times are estimated from the neighboring steps", file=sys.stderr)
    if n_outside:
        print(f"note: {n_outside} output step(s) are outside the steps in {dat} and were skipped",
              file=sys.stderr)


if __name__ == "__main__":
    main()

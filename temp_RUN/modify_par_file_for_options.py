import toml

import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--file_name", "-f", type=str, help="Where the file is")
parser.add_argument("--output_name", "-o", type=str, help="Where the file is")

# we want to vary grain size
parser.add_argument("--grain_size", "-g", type=int, help="Grain size to set")
parser.add_argument("--max_depth", "-m", type=int, help="Max depth to set")

args = parser.parse_args()

with open(args.file_name, "r") as f:
    data = toml.load(f)


if args.grain_size is not None:
    data["BSSN_DENDRO_GRAIN_SZ"] = args.grain_size

if args.max_depth is not None:
    data["BSSN_MAXDEPTH"] = args.max_depth


with open(args.output_name, "w") as f:
    f.write(toml.dumps(data))

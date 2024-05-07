import toml
import copy

with open("bssnq1_BASELINE.par.toml", "r") as f:
    base_options = toml.load(f)

print(base_options)

# need to change these values
blosc_options = [3, 6, 9]
zfp_options = [
    ("accuracy", 1e-5),
    ("accuracy", 1e-4),
    ("accuracy", 1e-3),
    ("accuracy", 1e-2),
    ("precision", 32),
    ("precision", 28),
    ("precision", 24),
    ("precision", 20),
    ("precision", 16),
]

BASH_PREFIX = "mpirun -np 12 ./bssnSolver ../../parstest/"

bash_script = ""

for send_type in ["double", "float"]:
    for compression_type in ["blosc", "zfp"]:
        temp_options = copy.deepcopy(base_options)

        temp_options["BSSN_COMPRESSION_SEND_TYPE"] = send_type

        if compression_type == "blosc":
            for blosc_option in blosc_options:
                temp_options["BSSN_COMPRESSION_OPTIONS"]["BLOSC_LEVEL"] = blosc_option
                temp_options["BSSN_COMPRESSION_MODE"] = "blosc"
                file_prefix = "bssnq1_BLOSC_" + send_type + "_Level" + str(blosc_option)
                temp_options["BSSN_VTU_FILE_PREFIX"] = file_prefix

                with open("parstest/" + file_prefix + ".par.toml", "w") as f:
                    toml.dump(temp_options, f)

                bash_script += (
                    BASH_PREFIX
                    + file_prefix
                    + ".par.toml | tee "
                    + file_prefix
                    + "_out.txt \n\n"
                )

        elif compression_type == "zfp":
            for ii, zfp_option in enumerate(zfp_options):
                temp_options["BSSN_COMPRESSION_MODE"] = "zfp"
                if zfp_option[0] == "accuracy":
                    temp_options["BSSN_COMPRESSION_OPTIONS"]["ZFP_ACCURACY"] = (
                        zfp_option[1]
                    )
                    temp_options["BSSN_COMPRESSION_OPTIONS"]["ZFP_MODE"] = "accuracy"
                elif zfp_option[0] == "precision":
                    temp_options["BSSN_COMPRESSION_OPTIONS"]["ZFP_PRECISION"] = (
                        zfp_option[1]
                    )
                    temp_options["BSSN_COMPRESSION_OPTIONS"]["ZFP_MODE"] = "precision"

                file_prefix = "bssnq1_ZFP_" + send_type + "_" + zfp_option[0] + str(ii)
                temp_options["BSSN_VTU_FILE_PREFIX"] = file_prefix

                with open("parstest/" + file_prefix + ".par.toml", "w") as f:
                    toml.dump(temp_options, f)

                bash_script += (
                    BASH_PREFIX
                    + file_prefix
                    + ".par.toml | tee "
                    + file_prefix
                    + "_out.txt \n\n"
                )


with open("runall.sh", "w") as f:
    bash_script = "#!/bin/bash\n\ntrap 'exit 120' INT\n\n" + bash_script
    f.write(bash_script)

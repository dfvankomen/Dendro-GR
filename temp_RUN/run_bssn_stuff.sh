# MAX_PROCS=(4 8 16 32 64 128 256 512 1024 2048 4096 8192)
MAX_PROCS=(4 8 16 32 64 128)
GRAIN_SIZES=(20 30 40 50 60 70 80 90 100)
MAX_DEPTHS=(13 14 15)
# MAX_DEPTHS=(10)

for maxdepth in "${MAX_DEPTHS[@]}"; do
  for procsize in "${MAX_PROCS[@]}"; do
    for gsz in "${GRAIN_SIZES[@]}"; do
      echo "Executing command with procsize=$procsize and grainsize=$gsz and maxdepth=$maxdepth"

      python modify_par_file_for_options.py --file_name=q1.par.toml --output_name=pars_to_use.toml --grain_size="$gsz" --max_depth="$maxdepth"

      output_file="output_np${procsize}_gsz${gsz}_maxdepth${maxdepth}.txt"
      #
      mpirun --oversubscribe -np "$procsize" ./bssnSolver pars_to_use.toml >"$output_file"
    done
  done
done

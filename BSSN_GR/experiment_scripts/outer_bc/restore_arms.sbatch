#!/bin/bash -l
#SBATCH --job-name=bssnobc
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.err
#SBATCH --partition=m12
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --mem=480G
#SBATCH --time=05:00:00
#
# Continue arms of an earlier run_arms job from their latest checkpoints to a new end time.
#   FROM=<job id> TEND=60.0 sbatch restore_arms.sbatch
#   FROM=<job id> TEND=60.0 ARMS="ghost_A6" sbatch --nodes=2 restore_arms.sbatch
# The job cancels itself unless each restarted arm reproduces a constraint row of the
# original run exactly (a wrong restart or a different binary would not).
SRC=$HOME/Dendro-GR
BUILD=$HOME/build_obc                    # must be the binary that wrote the checkpoints
ARMS=${ARMS:-"ghost_A6 ghost_E6"}
: "${FROM:?set FROM=<job id of the run to continue>}" "${TEND:?set TEND=<new end time>}"

source ~/gr_env.sh
ulimit -c 0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
SRCROOT=$HOME/nobackup/autodelete/bssnobc_$FROM
RUNROOT=$HOME/nobackup/autodelete/bssnobc_$SLURM_JOB_ID
BIN=$BUILD/BSSN_GR/bssnSolver
mkdir -p $RUNROOT
echo "host $(hostname) date $(date -u +%FT%TZ) from $FROM to t = $TEND, bin sha $(sha256sum $BIN | cut -c1-16)"

HOSTS=($(scontrol show hostnames $SLURM_JOB_NODELIST))
restore_arm() {  # arm host_pair
    local dir=$RUNROOT/$1
    rsync -a --exclude vtu $SRCROOT/$1/ $dir/ && mkdir -p $dir/vtu
    rm -f $dir/exit.txt && mv $dir/run.log $dir/run_before_restore.log
    sed -i -e "s/^BSSN_RESTORE_SOLVER = .*/BSSN_RESTORE_SOLVER = 1/" \
           -e "s/^BSSN_RK_TIME_END = .*/BSSN_RK_TIME_END = $TEND/" $dir/q1.par.toml
    (cd $dir && srun --mpi=pmix --exact -N2 -n256 -w $2 $BIN q1.par.toml > run.log 2>&1
     echo "exit=$?" > exit.txt)
}
i=0
for arm in $ARMS; do
    restore_arm $arm "${HOSTS[$i]},${HOSTS[$((i + 1))]}" &
    i=$((i + 2))
done

# Rows after the checkpoint are recomputed; the first one must match the original.
sleep 900
for arm in $ARMS; do
    norms=$RUNROOT/$arm/dat/dgr_ConstraintNorms.dat
    step=$(grep '^[0-9]' $SRCROOT/$arm/dat/dgr_ConstraintNorms.dat | tail -1 | cut -f1)
    copies=$(grep -c "^$step	" $norms); distinct=$(grep "^$step	" $norms | sort -u | wc -l)
    echo "$arm: step $step rows $copies distinct $distinct"
    [ "$copies" -ge 2 ] && [ "$distinct" -eq 1 ] \
        || { echo "ERROR: $arm restart does not reproduce step $step"; scancel $SLURM_JOB_ID; }
done
wait
for arm in $ARMS; do
    echo "$arm $(cat $RUNROOT/$arm/exit.txt) last row: $(tail -1 $RUNROOT/$arm/dat/dgr_Constraints.dat | cut -f1-3)"
done
echo "all done $(date -u +%T)"; du -sh $RUNROOT

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
# Outer-BC arms of the q1 BBH from fresh initial data, 2 nodes (256 ranks) per arm.
#   sbatch run_arms.sbatch                       # ARMS default: ghost_A6 ghost_E6
#   ARMS="ghost_A6 ghost_E6 perfield_A6" sbatch --nodes=6 run_arms.sbatch
# Each arm name is a parfile pars/<arm>.par.toml. Set --nodes to 2 x the number of arms.
# Edit the four paths below once for your account.
SRC=$HOME/Dendro-GR                     # this repository, branch exp/bssn-outer-bc
BUILD=$HOME/build_obc                    # build directory (made on first use)
DENDROLIB=$HOME/dendrolib                # dendrolib, branch experimental
ID=$HOME/id/tp_q001_tpid_sol.bin        # TwoPunctures q1 initial data
ARMS=${ARMS:-"ghost_A6 ghost_E6"}

source ~/gr_env.sh                       # modules: compilers, MPI, ...
ulimit -c 0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
KIT=$SRC/BSSN_GR/experiment_scripts/outer_bc
RUNROOT=$HOME/nobackup/autodelete/bssnobc_$SLURM_JOB_ID
mkdir -p $RUNROOT
echo "host $(hostname) date $(date -u +%FT%TZ) commit $(git -C $SRC rev-parse --short HEAD)"
[ -f $ID ] || { echo "ERROR: missing initial data $ID"; exit 2; }

if [ ! -x $BUILD/BSSN_GR/bssnSolver ]; then
    cmake -S $SRC -B $BUILD -DCMAKE_BUILD_TYPE=Release -DWITH_CUDA=OFF -DCPU_ARCH=znver3 \
        -DDENDRO_dendrolib_DIR=$DENDROLIB -DDENDRO_USE_NEW_DERIVS=ON > $BUILD.cfg.log 2>&1 \
        || { echo "CONFIGURE FAILED"; tail -30 $BUILD.cfg.log; exit 2; }
    cmake --build $BUILD --target bssnSolver -j 128 > $BUILD.bld.log 2>&1 \
        || { echo "BUILD FAILED"; grep -m5 -B2 -A8 error $BUILD.bld.log; exit 2; }
fi
BIN=$BUILD/BSSN_GR/bssnSolver
echo "bssnSolver sha $(sha256sum $BIN | cut -c1-16)"

HOSTS=($(scontrol show hostnames $SLURM_JOB_NODELIST))
run_arm() {  # arm host_pair
    local dir=$RUNROOT/$1
    mkdir -p $dir/cp $dir/bah $dir/dat $dir/vtu
    cp $KIT/pars/$1.par.toml $dir/q1.par.toml
    cp $ID $dir/
    (cd $dir && srun --mpi=pmix --exact -N2 -n256 -w $2 $BIN q1.par.toml > run.log 2>&1
     echo "exit=$?" > exit.txt)
}
i=0
for arm in $ARMS; do
    run_arm $arm "${HOSTS[$i]},${HOSTS[$((i + 1))]}" &
    i=$((i + 2))
done

# The startup line must name the configuration each arm was meant to run.
sleep 300
for arm in $ARMS; do
    echo "$arm: $(sed 's/\x1b\[[0-9;]*m//g' $RUNROOT/$arm/run.log | grep -a -m1 'Outer BC')"
done
wait
for arm in $ARMS; do
    echo "$arm $(cat $RUNROOT/$arm/exit.txt) last row: $(tail -1 $RUNROOT/$arm/dat/dgr_Constraints.dat | cut -f1-3)"
done
echo "all done $(date -u +%T)"; du -sh $RUNROOT

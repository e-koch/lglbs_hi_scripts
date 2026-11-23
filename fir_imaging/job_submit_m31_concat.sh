#!/bin/bash
#SBATCH --time=216:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=m31-concat-configs
#SBATCH --output=m31-concat-configs-%J.out
#SBATCH --mail-user=ekoch@ualberta.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Concatenate the per-config M31 staged MSs into a combined MS. One job per combination
# so they can run in parallel:
#
#   sbatch job_submit_m31_concat.sh ABCD    ->  m31_A+B+C+D_<line>.ms.tar
#   sbatch job_submit_m31_concat.sh CD      ->  m31_C+D_<line>.ms.tar
#
# The combined tar is written directly into staged_measurement_sets. By default an existing tar of
# the same name is replaced; to skip lines whose combined tar already exists instead:
#   M31_CONCAT_OVERWRITE=0 sbatch job_submit_m31_concat.sh CD
#
# Each job runs in its own scratch dir: /home/ekoch/scratch/VLAXL_imaging/m31_concat_<combo>_<jobid>

set -e

COMBO=$1
case "$COMBO" in
    ABCD|A+B+C+D) COMBO=ABCD ;;
    CD|C+D)       COMBO=CD ;;
    *) echo "Usage: sbatch job_submit_m31_concat.sh <ABCD|CD>"; exit 1 ;;
esac

WORKDIR=/home/ekoch/scratch/VLAXL_imaging/m31_concat_${COMBO}_${SLURM_JOB_ID}
CASA_EXECUTABLE="/home/ekoch/casa-6.7.6-14-py3.12.el9/bin/casa"
CASA_SCRIPT="/home/ekoch/lglbs_hi_scripts/fir_imaging/run_m31_concat_configs.py"

module load StdEnv
module load qt
source /home/ekoch/.bashrc
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/usr/
export CASALD_LIBRARY_PATH=$LD_LIBRARY_PATH

mkdir -p $WORKDIR
cd $WORKDIR
echo "Work dir: $WORKDIR"

# Read by run_m31_concat_configs.py (casa -c does not pass script args reliably).
export M31_CONCAT_COMBO=$COMBO
export M31_CONCAT_LOCAL_DIR=$WORKDIR

# CASA >= 6.7 no longer has --rcdir or reads an rc file from it. Memory limits are
# set from Python instead, via a per-job startup file passed with --startupfile.
# casalog.setMemoryTotal() takes MB; casalog.setMemoryFraction() takes an INTEGER percent.
# SLURM_MEM_PER_NODE is already in MB. ~/.casa/config.py is read automatically.
make_casa_startup () {
    local alloc_mb=${SLURM_MEM_PER_NODE:-$(( ${SLURM_MEM_PER_CPU:-4000} * ${SLURM_CPUS_PER_TASK:-1} ))}
    local mem_mb=$(( alloc_mb * 9 / 10 ))      # ~10% headroom for OS / xvfb
    CASA_CACHEDIR=${SLURM_TMPDIR:-/tmp}/casa_cache_${SLURM_JOB_ID}
    mkdir -p $CASA_CACHEDIR
    CASA_STARTUP=$CASA_CACHEDIR/startup.py
    cat > $CASA_STARTUP <<EOF
from casatasks import casalog
casalog.setMemoryTotal(${mem_mb})
casalog.setMemoryFraction(99)
EOF
    echo "CASA cachedir $CASA_CACHEDIR: memory ${mem_mb} MB of ${alloc_mb} MB allocated, memfrac 99%"
}
make_casa_startup

xvfb-run -a $CASA_EXECUTABLE --cachedir $CASA_CACHEDIR --startupfile $CASA_STARTUP \
    --nologger --nogui --log2term -c $CASA_SCRIPT

seff $SLURM_JOB_ID || true

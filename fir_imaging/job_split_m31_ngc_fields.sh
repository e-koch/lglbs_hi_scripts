#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=m31-split-ngc
#SBATCH --output=m31-split-ngc-%J.out
#SBATCH --mail-user=ekoch@ualberta.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Split NGC185 / NGC205 out of the staged per-config M31 tars (see run_m31_split_ngc_fields.py).
# Needed for D tars staged before the staging key was changed to select only M31* fields; run it
# BEFORE job_submit_m31_concat.sh so the concats are built from clean per-config tars.
#
#   sbatch job_split_m31_ngc_fields.sh                   # D config, all lines (contsub + nocontsub)
#
# Optional env vars (inherited by sbatch):
#   M31_SPLIT_CONFIGS=A,B,C,D     configs to process (default D)
#   M31_SPLIT_LINES=hilores,...   line names, exact match (default all)
#   M31_SPLIT_DRY_RUN=1           only read each tar's FIELD table and report; writes nothing.
#                                 e.g. M31_SPLIT_DRY_RUN=1 M31_SPLIT_CONFIGS=A,B,C,D \
#                                      sbatch --time=02:00:00 job_split_m31_ngc_fields.sh
#   M31_SPLIT_STAGED_DIR=...      point at a test copy of staged_measurement_sets
#
# Output: ngc185_D_<line>[_nocontsub].ms.tar and ngc205_D_<...>.ms.tar in staged_measurement_sets,
# and m31_D_<line>[_nocontsub].ms.tar replaced by the M31*-only version. A tar that is already
# clean is left alone, so the job is safe to repeat. Peak scratch is ~2.2x the largest D MS
# (original + the split copies). --time is a guess (16 D tars, one at a time): check the log.
#
# Each job runs in its own scratch dir: /home/ekoch/scratch/VLAXL_imaging/m31_ngc_split_<jobid>

set -e

WORKDIR=/home/ekoch/scratch/VLAXL_imaging/m31_ngc_split_${SLURM_JOB_ID}
CASA_EXECUTABLE="/home/ekoch/casa-6.7.6-14-py3.12.el9/bin/casa"
CASA_SCRIPT="/home/ekoch/lglbs_hi_scripts/fir_imaging/run_m31_split_ngc_fields.py"

module load StdEnv
module load qt
source /home/ekoch/.bashrc
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/usr/
export CASALD_LIBRARY_PATH=$LD_LIBRARY_PATH

mkdir -p $WORKDIR
cd $WORKDIR
echo "Work dir: $WORKDIR"

# Read by run_m31_split_ngc_fields.py (casa -c does not pass script args reliably).
export M31_SPLIT_LOCAL_DIR=$WORKDIR

# Memory limits via a per-job startup file (same setup as job_submit_m31_concat.sh).
make_casa_startup () {
    local alloc_mb=${SLURM_MEM_PER_NODE:-$(( ${SLURM_MEM_PER_CPU:-4000} * ${SLURM_CPUS_PER_TASK:-1} ))}
    local mem_mb=$(( alloc_mb * 9 / 10 ))
    CASA_CACHEDIR=${SLURM_TMPDIR:-/tmp}/casa_cache_${SLURM_JOB_ID}
    mkdir -p $CASA_CACHEDIR
    CASA_STARTUP=$CASA_CACHEDIR/startup.py
    cat > $CASA_STARTUP <<EOT
from casatasks import casalog
casalog.setMemoryTotal(${mem_mb})
casalog.setMemoryFraction(99)
EOT
    echo "CASA cachedir $CASA_CACHEDIR: memory ${mem_mb} MB of ${alloc_mb} MB allocated, memfrac 99%"
}
make_casa_startup

xvfb-run -a $CASA_EXECUTABLE --cachedir $CASA_CACHEDIR --startupfile $CASA_STARTUP \
    --nologger --nogui --log2term -c $CASA_SCRIPT

seff $SLURM_JOB_ID || true

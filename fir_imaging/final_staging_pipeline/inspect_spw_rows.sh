#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --job-name=inspect_spw_rows
#SBATCH --output=inspect_spw_rows-%J.out
#SBATCH --mail-user=ekoch@ualberta.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Row / flag statistics per (observation, SPW) of a staged MS tar (see inspect_spw_rows.py).
# The whole tar is extracted into a per-job scratch dir and removed afterwards, so allow scratch for
# one tar (oh1665 ABCD is ~290 GB). Run on the ORIGINAL tar, before job 2 replaces it.
#
#   sbatch inspect_spw_rows.sh                                  # m33_A+B+C+D_oh1665
#   INSPECT_TARS="m33_A+B+C+D_oh1667" sbatch inspect_spw_rows.sh
#   INSPECT_TARS="m33_A+B+C+D_oh1665,m33_A+B+C+D_oh1667,m33_C+D_oh1667" sbatch --time=24:00:00 inspect_spw_rows.sh
#
# --time is a guess (extract + one pass over FLAG/WEIGHT per tar): check the log of the first run.
# Other env vars: INSPECT_STAGED_DIR, INSPECT_OUTDIR (CSV per tar; default: the submit dir),
# INSPECT_LINE_MHZ, INSPECT_HALFWIDTH_KMS, INSPECT_KEEP_MS.

set -e
CASA_EXECUTABLE="/home/ekoch/casa-6.7.6-14-py3.12.el9/bin/casa"
SCRIPT="/home/ekoch/lglbs_hi_scripts/fir_imaging/final_staging_pipeline/inspect_spw_rows.py"

module load StdEnv
module load qt
source /home/ekoch/.bashrc
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/usr/
export CASALD_LIBRARY_PATH=$LD_LIBRARY_PATH

export INSPECT_WORKDIR=/home/ekoch/scratch/VLAXL_imaging/inspect_spw_rows_${SLURM_JOB_ID}
export INSPECT_OUTDIR=${INSPECT_OUTDIR:-$SLURM_SUBMIT_DIR}
mkdir -p $INSPECT_WORKDIR
cd $INSPECT_WORKDIR
echo "Work dir: $INSPECT_WORKDIR"

xvfb-run -a $CASA_EXECUTABLE --nologger --nogui --log2term -c $SCRIPT

seff $SLURM_JOB_ID || true

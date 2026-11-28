#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=inspect_spws
#SBATCH --output=inspect_spws-%J.out
#SBATCH --mail-user=ekoch@ualberta.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# List the SPWs of staged MS tars without extracting the data (see inspect_spws_from_tar.py).
# Run on the ORIGINAL tars in staged_measurement_sets, before job 2 replaces them.
#
#   sbatch inspect_spws_from_tar.sh                       # the M33 OH tars (A+B+C+D and C+D)
#   INSPECT_TARS="m33_A+B+C+D_oh1665,m33_A+B+C+D_oh1667" sbatch inspect_spws_from_tar.sh
#
# Other env vars: INSPECT_STAGED_DIR, INSPECT_OUTDIR (CSV per tar; default: the submit dir),
# INSPECT_LINE_MHZ, INSPECT_HALFWIDTH_KMS. Reading only the SPW subtables still scans the tar's
# headers, so allow a few minutes per tar on project space.

set -e
CASA_EXECUTABLE="/home/ekoch/casa-6.7.6-14-py3.12.el9/bin/casa"
SCRIPT="/home/ekoch/lglbs_hi_scripts/fir_imaging/final_staging_pipeline/inspect_spws_from_tar.py"

module load StdEnv
module load qt
source /home/ekoch/.bashrc
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/usr/
export CASALD_LIBRARY_PATH=$LD_LIBRARY_PATH
export INSPECT_OUTDIR=${INSPECT_OUTDIR:-$SLURM_SUBMIT_DIR}

xvfb-run -a $CASA_EXECUTABLE --nologger --nogui --log2term -c $SCRIPT

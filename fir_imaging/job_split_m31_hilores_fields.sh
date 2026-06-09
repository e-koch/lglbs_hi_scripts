#!/bin/bash
#SBATCH --job-name=split_m31_hilores
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=split_m31_hilores-%j.out
#SBATCH --mail-user=ekoch@ualberta.ca
#SBATCH --mail-type=END,FAIL

set -e

cd /home/ekoch/scratch/VLAXL_imaging/imaging/m31

CASA_EXECUTABLE="/home/ekoch/casa-6.6.1-17-pipeline-2024.1.0.8/bin/casa"
SCRIPT="/home/ekoch/lglbs_hi_scripts/fir_imaging/split_m31_hilores_fields.py"

module load StdEnv
module load qt

# Random sleep to avoid job start conflicts
xvfb-run -a $CASA_EXECUTABLE --rcdir ~/.casa --nologger --nogui --log2term -c $SCRIPT

echo "M31 hilores field split complete"

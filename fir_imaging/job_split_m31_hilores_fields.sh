#!/bin/bash
#SBATCH --job-name=split_m31_tiles
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=0-3
#SBATCH --output=split_m31_tiles-%A_%a.out
#SBATCH --mail-user=ekoch@ualberta.ca
#SBATCH --mail-type=END,FAIL

# Split the staged M31 MS into the 4 strip-tile MSes (one array task per
# tile; see m31_tile_definitions.py). Runs after staging, before per-tile
# imaging. Optional args: $1 = config (default A+B+C+D),
#                         $2 = product (default hilores)
#
# Submit:  sbatch ~/lglbs_hi_scripts/fir_imaging/job_split_m31_hilores_fields.sh
#     or:  sbatch ... job_split_m31_hilores_fields.sh A+B+C+D himidres

set -e

TILES=(m31_1 m31_2 m31_3 m31_4)
TILE=${TILES[$SLURM_ARRAY_TASK_ID]}
CONFIG=${1:-A+B+C+D}
PRODUCT=${2:-hilores}

module load StdEnv
module load qt

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

CASA_EXECUTABLE="/home/ekoch/casa-6.6.1-17-pipeline-2024.1.0.8/bin/casa"
SCRIPT="/home/ekoch/lglbs_hi_scripts/fir_imaging/split_m31_hilores_fields.py"

cd /home/ekoch/scratch/VLAXL_imaging/imaging/m31

# Stagger array-task starts so 4 CASA instances don't collide on startup
python3 -c "import time, random; time.sleep(random.randint(2, 60))"

xvfb-run -a $CASA_EXECUTABLE --rcdir ~/.casa --nologger --nogui --log2term \
    -c $SCRIPT $TILE $CONFIG $PRODUCT

echo "M31 tile split complete: $TILE ($CONFIG $PRODUCT)"
seff $SLURM_JOB_ID || true

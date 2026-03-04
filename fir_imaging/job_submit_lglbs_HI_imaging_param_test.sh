#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=lglbs-HI-test-imaging-params-%J
#SBATCH --output=lglbs-HI-test-imaging-params-%J.out
#SBATCH --mail-user=ekoch@ualberta.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# For now, just pass the target name from the cmd line with sbatch.
# Will want the ability to pass an array of target names later.
export target_name=$1
export config_name=$2
export product_name=$3

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load StdEnv
module load qt

source /home/ekoch/.bashrc

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/usr/

export CASALD_LIBRARY_PATH=$LD_LIBRARY_PATH


export casa_executable="/home/ekoch/casa-6.6.1-17-pipeline-2024.1.0.8/bin/casa"
export casa_script="/home/ekoch/lglbs_hi_scripts/fir_imaging/run_lglbs_HI_test_imaging_params.py"


# Step 2: Run CASA to do staging.
xvfb-run -a $casa_executable --rcdir ~/.casa --nologger --nogui --log2term -c $casa_script $target_name $config_name $product_name


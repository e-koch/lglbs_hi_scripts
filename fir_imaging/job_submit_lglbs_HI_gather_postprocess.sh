#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --job-name=lglbs-hi-gather-postprocess-%J
#SBATCH --output=lglbs-hi-gather-postprocess-%J.out
#SBATCH --mail-user=ekoch@ualberta.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# For now, just pass the target name from the cmd line with sbatch.
# Will want the ability to pass an array of target names later.

export this_galaxy=$1
# export this_line_product='hilores'
export this_line_product='himidres'
export this_config='A+B+C+D'
export this_chunksize=20
export this_idx=${SLURM_ARRAY_TASK_ID}


echo $SLURM_CPUS_PER_TASK

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load StdEnv
module load qt

source /home/ekoch/.bashrc

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/usr/

export CASALD_LIBRARY_PATH=$LD_LIBRARY_PATH


# Ensure no time overlap in job start times
python3 -c "import time, random; time.sleep(random.randint(2, 120))"

export casa_executable="/home/ekoch/casa-6.6.1-17-pipeline-2024.1.0.8/bin/casa"
export casa_script="/home/ekoch/lglbs_hi_scripts/fir_imaging/run_lglbs_HI_gather_postprocess.py"
# export casa_job_config_file=/home/ekoch/lglbs_hi_scripts/fir_imaging/line_staging_imaging.${SLURM_ARRAY_TASK_ID}.jobconfig.txt


export script_args="$this_galaxy $this_config $this_line_product $this_chunksize"
echo "Args passed to script: $script_args"
xvfb-run -a $casa_executable --rcdir ~/.casa --nologger --nogui --log2term -c $casa_script $script_args

#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=16
#SBATCH --array=0-1
#SBATCH --job-name=lglbs-staging-%A-%a
#SBATCH --output=lglbs-staging-%A-%a.out
#SBATCH --mail-user=ekoch@ualberta.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# For now, just pass the target name from the cmd line with sbatch.
# Will want the ability to pass an array of target names later.

export this_galaxy=$1
export this_line_product='hilores'
export this_config='A+B+C+D'
export this_chunksize=25
export this_idx=${SLURM_ARRAY_TASK_ID}

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load StdEnv
module load qt

source /home/ekoch/.bashrc

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/usr/

export CASALD_LIBRARY_PATH=$LD_LIBRARY_PATH


# Ensure no time overlap in job start times
python3 -c "import time, random; time.sleep(random.randint(2, 120))"

export data_path="/home/ekoch/scratch/VLAXL_imaging/MeasurementSets/"

export casa_executable="/home/ekoch/casa-6.6.1-17-pipeline-2024.1.0.8/bin/casa"
export casa_script="/home/ekoch/lglbs_hi_scripts/fir_imaging/run_lglbs_HI_imaging.py"
# export casa_job_config_file=/home/ekoch/lglbs_hi_scripts/fir_imaging/line_staging_imaging.${SLURM_ARRAY_TASK_ID}.jobconfig.txt


# Step 1: Check if we need to untar the MSs
# Set variable to current working directory
curr_dir=$(pwd)

cd $data_path

if [ ! -d "$target_name" ]; then
    echo "No folder $target_name. Restore from project/nearline first with globus"
fi

cd $target_name
# Loop through all tar files, untar them, and delete them.

for filename in *.tar; do
    echo "Untarring $filename"
    tar -xf $filename
    rm $filename
done

cd $curr_dir


# Step 2: Run CASA to do staging.
export script_args="$this_galaxy $this_line_product $this_config $this_chunksize $this_idx"
echo "Args passed to script: $script_args"
xvfb-run -a $casa_executable --rcdir ~/.casa --nologger --nogui --log2term -c $casa_script $script_args


#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --job-name=lglbs-products-%J
#SBATCH --output=lglbs-products-%J.out
#SBATCH --mail-user=ekoch@ualberta.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


source /home/ekoch/.bashrc


# python 3.11 install
module load StdEnv/2023

# ipython/numpy/scipy/mpl
module load scipy-stack

# Setup python environment environment:
source $HOME/phangs/bin/activate


# For now, just pass the target name from the cmd line with sbatch.
# Will want the ability to pass an array of target names later.

export this_galaxy=$1
# export this_line_product='hilores'
export this_line_product='himidres'
export this_config='A+B+C+D'

export python_script=$HOME/lglbs_hi_scripts/fir_imaging/run_lglbs_HI_products.py

export script_args="$this_galaxy $this_config $this_line_product"
echo "Args passed to script: $script_args"

echo "python $python_script $script_args"

python $python_script $script_args

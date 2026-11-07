#!/bin/bash
#SBATCH --job-name=test_img_single
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=test_img_single-%J.out
#SBATCH --mail-user=ekoch@ualberta.ca
#SBATCH --mail-type=END,FAIL

# Re-run ONLY the test-imaging step of job 1 (run_test_imaging.py) on one specific,
# already-extracted MS. Nothing else from job1 (extract, combine_spws, weight-outlier
# flagging, uvcontsub) is run, so the MS must already be in the state you want to image.
#
# Submit:  sbatch test_imaging_single_ms.sh <ms_dir> <ms_name> [outdir]
#     e.g. sbatch ~/lglbs_hi_scripts/fir_imaging/final_staging_pipeline/test_imaging_single_ms.sh \
#              /home/ekoch/scratch/VLAXL_imaging/ic1613ctr_12345678 ic1613ctr_OH1667_A.ms
#
# <ms_dir>   directory containing the MS
# <ms_name>  MS directory name (with or without trailing slash)
# [outdir]   where test_<field>_<config>.{psf,residual} and test_imaging_summary.txt go.
#            Default: <ms_dir>/<ms_name minus .ms>_products/test_images (the job1 location,
#            so summarize_stage1_products.py picks the results up; existing images with the
#            same names are overwritten).

set -e

MS_DIR=${1%/}
MS_NAME=$(basename "${2%/}")
if [ -z "$MS_DIR" ] || [ -z "$MS_NAME" ]; then
    echo "Usage: sbatch test_imaging_single_ms.sh <ms_dir> <ms_name> [outdir]"; exit 1
fi

MS=$MS_DIR/$MS_NAME
if [ ! -d "$MS" ]; then
    echo "MS not found: $MS"; exit 1
fi

BASE=${MS_NAME%.ms}
OUTDIR=${3:-$MS_DIR/${BASE}_products/test_images}

PIPE=/home/ekoch/lglbs_hi_scripts/fir_imaging/final_staging_pipeline
CASA_EXECUTABLE="/home/ekoch/casa-6.7.6-14-py3.12.el9/bin/casa"

module load StdEnv
module load qt
source /home/ekoch/.bashrc
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/usr/
export CASALD_LIBRARY_PATH=$LD_LIBRARY_PATH

mkdir -p "$OUTDIR"
cd "$OUTDIR"
echo "MS:      $MS"
echo "Outdir:  $OUTDIR"

# Same CASA memory setup as job1_prep_and_test.sh (CASA >= 6.7: set limits from a
# per-job startup file; setMemoryTotal takes MB, setMemoryFraction an integer percent).
alloc_mb=${SLURM_MEM_PER_NODE:-$(( ${SLURM_MEM_PER_CPU:-4000} * ${SLURM_CPUS_PER_TASK:-1} ))}
mem_mb=$(( alloc_mb * 9 / 10 ))
CASA_CACHEDIR=${SLURM_TMPDIR:-/tmp}/casa_cache_${SLURM_JOB_ID}
mkdir -p "$CASA_CACHEDIR"
CASA_STARTUP=$CASA_CACHEDIR/startup.py
cat > "$CASA_STARTUP" <<EOF
from casatasks import casalog
casalog.setMemoryTotal(${mem_mb})
casalog.setMemoryFraction(99)
EOF
echo "CASA cachedir $CASA_CACHEDIR: memory ${mem_mb} MB of ${alloc_mb} MB allocated, memfrac 99%"

echo "$(date): test imaging $MS_NAME"
xvfb-run -a $CASA_EXECUTABLE --cachedir "$CASA_CACHEDIR" --startupfile "$CASA_STARTUP" \
    --nologger --nogui --log2term -c $PIPE/run_test_imaging.py "$MS" "$OUTDIR"

echo "$(date): done. Results in $OUTDIR"
seff $SLURM_JOB_ID || true

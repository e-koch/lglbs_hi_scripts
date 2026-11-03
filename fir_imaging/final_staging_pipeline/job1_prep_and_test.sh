#!/bin/bash
#SBATCH --job-name=stage_qa_job1
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=stage_qa_job1-%J.out
#SBATCH --mail-user=ekoch@ualberta.ca
#SBATCH --mail-type=END,FAIL

# Job 1 of the final staging pipeline. For one target, for each staged MS tar:
#   extract -> find_weight_outliers (+apply flags, flag summary)
#           -> uvcontsub (skipped for *nocontsub*) -> single-channel test imaging
# Everything lands in /home/ekoch/scratch/VLAXL_imaging/<target>_<jobid>/ so the
# outputs can be inspected before job2_archive_and_cleanup.sh is submitted.
#
# Submit:  sbatch job1_prep_and_test.sh <target>
#     e.g. sbatch ~/lglbs_hi_scripts/fir_imaging/final_staging_pipeline/job1_prep_and_test.sh ic1613ctr

set -e

TARGET=$1
if [ -z "$TARGET" ]; then
    echo "Usage: sbatch job1_prep_and_test.sh <target>"; exit 1
fi

STAGED_DIR=${STAGED_DIR:-/home/ekoch/projects/rrg-eros-ab/ekoch/VLAXL/staged_measurement_sets}
WORKDIR=/home/ekoch/scratch/VLAXL_imaging/${TARGET}_${SLURM_JOB_ID}
PIPE=/home/ekoch/lglbs_hi_scripts/fir_imaging/final_staging_pipeline
CASA_EXECUTABLE="/home/ekoch/casa-6.6.1-17-pipeline-2024.1.0.8/bin/casa"

module load StdEnv
module load qt
source /home/ekoch/.bashrc
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/usr/
export CASALD_LIBRARY_PATH=$LD_LIBRARY_PATH

mkdir -p $WORKDIR
cd $WORKDIR
echo "Work dir: $WORKDIR"
echo $WORKDIR > /home/ekoch/scratch/VLAXL_imaging/last_workdir_${TARGET}.txt

shopt -s nullglob
TARS=($STAGED_DIR/${TARGET}_*.ms.tar)
if [ ${#TARS[@]} -eq 0 ]; then
    echo "No staged tars matching $STAGED_DIR/${TARGET}_*.ms.tar"; exit 1
fi
echo "Found ${#TARS[@]} staged MS tar(s)"

# Build a per-job CASA rcdir whose memory settings match the SLURM allocation.
# system.resources.memory is in MB; SLURM_MEM_PER_NODE is already in MB.
make_casa_rcdir () {
    local mem_mb=${SLURM_MEM_PER_NODE:-$(( ${SLURM_MEM_PER_CPU:-4000} * ${SLURM_CPUS_PER_TASK:-1} ))}
    CASA_RCDIR=${SLURM_TMPDIR:-/tmp}/casa_rc_${SLURM_JOB_ID}
    mkdir -p $CASA_RCDIR
    [ -f ~/.casa/config.py ] && cp ~/.casa/config.py $CASA_RCDIR/
    printf 'system.resources.memory: %s\nsystem.resources.memfrac: 0.99\n' $mem_mb > $CASA_RCDIR/rc
    echo "CASA rcdir $CASA_RCDIR: memory ${mem_mb} MB, memfrac 0.99"
}
make_casa_rcdir

run_casa () {
    xvfb-run -a $CASA_EXECUTABLE --rcdir $CASA_RCDIR --nologger --nogui --log2term -c "$@"
}

for TAR in "${TARS[@]}"; do
    MS_NAME=$(basename $TAR .tar)          # <name>.ms
    BASE=${MS_NAME%.ms}
    PRODUCTS=$WORKDIR/${BASE}_products
    echo "=================================================="
    echo "$(date): $MS_NAME"
    echo "=================================================="

    # Extract straight from project space (strip leading path before *.ms)
    python3 -c "
import sys; sys.path.insert(0, '$PIPE')
from pipeline_common import extract_ms_tar
print(extract_ms_tar('$TAR', '$WORKDIR'))"
    mkdir -p $PRODUCTS

    run_casa $PIPE/find_weight_outliers.py $WORKDIR/$MS_NAME $PRODUCTS

    if [[ "$MS_NAME" == *nocontsub* ]]; then
        echo "nocontsub MS: skipping uvcontsub"
    else
        run_casa $PIPE/run_uvcontsub.py $WORKDIR/$MS_NAME
    fi

    run_casa $PIPE/run_test_imaging.py $WORKDIR/$MS_NAME $PRODUCTS/test_images
done

echo "Job 1 complete. Inspect $WORKDIR/*_products, then run:"
echo "  sbatch $PIPE/job2_archive_and_cleanup.sh $TARGET $WORKDIR"
seff $SLURM_JOB_ID || true

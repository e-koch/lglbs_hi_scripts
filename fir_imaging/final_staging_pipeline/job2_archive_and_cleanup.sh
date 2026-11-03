#!/bin/bash
#SBATCH --job-name=stage_qa_job2
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=stage_qa_job2-%J.out
#SBATCH --mail-user=ekoch@ualberta.ca
#SBATCH --mail-type=END,FAIL

# Job 2 of the final staging pipeline. Run AFTER checking job 1 outputs.
# For each updated MS in the Job 1 work dir: tar it, verify the tar, replace the
# tar in the project-space staged_measurement_sets, then remove the local MS.
# The *_products QA folders are kept unless PURGE=1.
#
# Submit:  sbatch job2_archive_and_cleanup.sh <target> [workdir]
#   workdir defaults to the one in scratch/VLAXL_imaging/last_workdir_<target>.txt
#   PURGE=1 sbatch --export=ALL,PURGE=1 ...  also deletes the work dir when done.

set -e

TARGET=$1
WORKDIR=$2
if [ -z "$TARGET" ]; then
    echo "Usage: sbatch job2_archive_and_cleanup.sh <target> [workdir]"; exit 1
fi
if [ -z "$WORKDIR" ]; then
    WORKDIR=$(cat /home/ekoch/scratch/VLAXL_imaging/last_workdir_${TARGET}.txt)
fi
[ -d "$WORKDIR" ] || { echo "Work dir not found: $WORKDIR"; exit 1; }

STAGED_DIR=${STAGED_DIR:-/home/ekoch/projects/rrg-eros-ab/ekoch/VLAXL/staged_measurement_sets}
PIPE=/home/ekoch/lglbs_hi_scripts/fir_imaging/final_staging_pipeline

module load StdEnv
cd $WORKDIR
echo "Work dir: $WORKDIR"

shopt -s nullglob
MSS=($WORKDIR/${TARGET}_*.ms)
[ ${#MSS[@]} -gt 0 ] || { echo "No MSs in $WORKDIR"; exit 1; }

FAILED=0
for MS in "${MSS[@]}"; do
    MS_NAME=$(basename $MS)
    TAR_LOCAL=$WORKDIR/${MS_NAME}.tar
    TAR_DEST=$STAGED_DIR/${MS_NAME}.tar
    echo "$(date): archiving $MS_NAME"

    if python3 -c "
import sys; sys.path.insert(0, '$PIPE')
from pipeline_common import make_ms_tar, verify_tar
make_ms_tar('$MS', '$TAR_LOCAL')
print('members:', verify_tar('$TAR_LOCAL', '$MS_NAME'))
" && cp $TAR_LOCAL ${TAR_DEST}.tmp && cmp -s $TAR_LOCAL ${TAR_DEST}.tmp \
      && mv -f ${TAR_DEST}.tmp $TAR_DEST; then
        rm -rf $MS $TAR_LOCAL
        echo "  replaced $TAR_DEST"
    else
        echo "  FAILED for $MS_NAME; leaving MS and original tar in place"
        rm -f ${TAR_DEST}.tmp
        FAILED=1
    fi
done

if [ $FAILED -ne 0 ]; then
    echo "Some MSs failed; work dir NOT cleaned: $WORKDIR"; exit 1
fi

if [ "${PURGE:-0}" = "1" ]; then
    cd /home/ekoch/scratch/VLAXL_imaging
    rm -rf $WORKDIR
    echo "Purged $WORKDIR"
else
    echo "Kept QA products in $WORKDIR (set PURGE=1 to delete)"
fi
seff $SLURM_JOB_ID || true

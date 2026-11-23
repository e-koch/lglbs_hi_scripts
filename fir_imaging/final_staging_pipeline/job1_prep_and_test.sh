#!/bin/bash
#SBATCH --job-name=stage_qa_job1
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --output=stage_qa_job1-%J.out
#SBATCH --mail-user=ekoch@ualberta.ca
#SBATCH --mail-type=END,FAIL

# Job 1 of the final staging pipeline. For one target, for each staged MS tar:
#   extract -> combine_spws (only if >1 SPW in use, e.g. OH lines with B0D0+A0C0)
#           -> find_weight_outliers (+apply flags, flag summary)
#           -> plot_flag_vs_channel (flagged % vs channel + frequency, PNG and CSV)
#           -> uvcontsub (skipped for *nocontsub*) -> single-channel test imaging
# then, once every MS is done, summarize_stage1_products (markdown summary + PNG gallery
# per MS and an overview table, in <workdir>/stage1_qa_reports).
# Everything lands in /home/ekoch/scratch/VLAXL_imaging/<target>_<jobid>/ so the
# outputs can be inspected before job2_archive_and_cleanup.sh is submitted.
#
# Submit:  sbatch job1_prep_and_test.sh <target>
#     e.g. sbatch ~/lglbs_hi_scripts/fir_imaging/final_staging_pipeline/job1_prep_and_test.sh ic1613ctr
#   MAKE_SUMMARY=0 sbatch job1_prep_and_test.sh <target>   skips the summary report at the end.
#   A failure in the flag plot or the summary report is only a warning; it does not stop the job.
#
# Default: every staged tar for <target>. For big targets (M31) select a subset with env vars
# (all optional, combinable; the tar names are <target>_<config>_<line>[_nocontsub].ms.tar):
#   CONFIGS=C+D                    comma list of array configs: A B C D C+D A+B+C+D
#   LINES=hilores,himidres         comma list of line names, exact match (hi21cm != hi21cm_0p8kms)
#   CONTSUB=yes|no                 yes = skip *_nocontsub tars, no = only *_nocontsub (default: both)
#   TAG=name                       workdir label (default: derived from the filters above, e.g. ABCD_hilores)
#   DRY_RUN=1                      list the selected tars + sizes + estimated scratch need, then exit.
#                                  Runs fine with plain `bash` on a login node.
#   e.g.  CONFIGS=A+B+C+D LINES=oh1612 sbatch job1_prep_and_test.sh m31
# Each submission gets its own workdir and needs its own job 2 (pass the workdir explicitly).
#
# Test imaging (one tclean per field x config) is sharded over TEST_SHARDS CASA processes
# (default 8, each in its own directory + CASA cache dir; TEST_SHARDS=1 = serial as before).
# Each shard gets 1/TEST_SHARDS of the CASA memory budget and CPUS/TEST_SHARDS OpenMP threads.

set -e

TARGET=$1
if [ -z "$TARGET" ]; then
    echo "Usage: [CONFIGS=..] [LINES=..] [CONTSUB=..] [DRY_RUN=1] sbatch job1_prep_and_test.sh <target>"; exit 1
fi

STAGED_DIR=${STAGED_DIR:-/home/ekoch/projects/rrg-eros-ab/ekoch/VLAXL/staged_measurement_sets}
PIPE=/home/ekoch/lglbs_hi_scripts/fir_imaging/final_staging_pipeline
CASA_EXECUTABLE="/home/ekoch/casa-6.7.6-14-py3.12.el9/bin/casa"

# ---- Tar selection ------------------------------------------------------------------------
shopt -s nullglob
ALL_TARS=($STAGED_DIR/${TARGET}_*.ms.tar)

select_tars () {
    local t b cfg line ncs
    for t in "${ALL_TARS[@]}"; do
        b=$(basename $t .ms.tar); b=${b#${TARGET}_}      # e.g. C+D_hi21cm_0p8kms_nocontsub
        cfg=${b%%_*}; line=${b#*_}; ncs=0                # configs never contain '_'
        if [[ $line == *_nocontsub ]]; then ncs=1; line=${line%_nocontsub}; fi
        [ -n "$CONFIGS" ] && [[ ",$CONFIGS," != *",$cfg,"* ]] && continue
        [ -n "$LINES" ]   && [[ ",$LINES,"   != *",$line,"* ]] && continue
        [ "${CONTSUB:-both}" = yes ] && [ $ncs = 1 ] && continue
        [ "${CONTSUB:-both}" = no ]  && [ $ncs = 0 ] && continue
        echo $t
    done
}
mapfile -t TARS < <(select_tars)

if [ ${#TARS[@]} -eq 0 ]; then
    echo "No staged tars matching $STAGED_DIR/${TARGET}_*.ms.tar" \
         "(CONFIGS='$CONFIGS' LINES='$LINES' CONTSUB='${CONTSUB:-both}')"; exit 1
fi
echo "Selected ${#TARS[@]} of ${#ALL_TARS[@]} staged MS tar(s)" \
     "(CONFIGS='${CONFIGS:-all}' LINES='${LINES:-all}' CONTSUB='${CONTSUB:-both}'):"
TOTAL_B=0; LARGEST_B=0
for t in "${TARS[@]}"; do
    sz=$(stat -c %s "$t" 2>/dev/null || stat -f %z "$t")
    TOTAL_B=$(( TOTAL_B + sz )); [ $sz -gt $LARGEST_B ] && LARGEST_B=$sz
    printf '  %6d G  %s\n' $(( sz / 1073741824 )) $(basename $t)
done
echo "Total $(( TOTAL_B / 1073741824 )) G; largest $(( LARGEST_B / 1073741824 )) G;" \
     "est. peak scratch ~$(( (TOTAL_B + 2 * LARGEST_B) / 1073741824 )) G" \
     "(all MSs stay in the workdir until job 2; uvcontsub needs +1x the current MS, combine_spws up to +2x)"
[ "${DRY_RUN:-0}" = "1" ] && exit 0

# Workdir label so concurrent / split submissions do not collide
if [ -z "${TAG+x}" ]; then
    TAG_PARTS=()
    if [ -n "$CONFIGS" ]; then c=${CONFIGS//+/}; TAG_PARTS+=("${c//,/-}"); fi
    if [ -n "$LINES" ]; then TAG_PARTS+=("${LINES//,/-}"); fi
    if [ -n "$CONTSUB" ] && [ "$CONTSUB" != both ]; then TAG_PARTS+=("contsub-$CONTSUB"); fi
    TAG=$(IFS=_; echo "${TAG_PARTS[*]}")
fi
WORKDIR=/home/ekoch/scratch/VLAXL_imaging/${TARGET}${TAG:+_$TAG}_${SLURM_JOB_ID}

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
# With a tag, also keep a per-tag pointer (the untagged file is overwritten by whichever job started last)
[ -n "$TAG" ] && echo $WORKDIR > /home/ekoch/scratch/VLAXL_imaging/last_workdir_${TARGET}_${TAG}.txt

# CASA >= 6.7 no longer has --rcdir or reads an rc file from it. Memory limits are
# set from Python instead, via a per-job startup file passed with --startupfile.
# casalog.setMemoryTotal() takes MB (checked: setMemoryTotal(1000) -> getMemoryTotal()
# = 1024000 kB); casalog.setMemoryFraction() takes an INTEGER percent.
# SLURM_MEM_PER_NODE is already in MB. ~/.casa/config.py is read automatically.
make_casa_startup () {
    ALLOC_MB=${SLURM_MEM_PER_NODE:-$(( ${SLURM_MEM_PER_CPU:-4000} * ${SLURM_CPUS_PER_TASK:-1} ))}
    local alloc_mb=$ALLOC_MB
    local mem_mb=$(( alloc_mb * 9 / 10 ))      # ~10% headroom for OS / xvfb / non-imaging use
    CASA_CACHEDIR=${SLURM_TMPDIR:-/tmp}/casa_cache_${SLURM_JOB_ID}
    mkdir -p $CASA_CACHEDIR
    CASA_STARTUP=$CASA_CACHEDIR/startup.py
    cat > $CASA_STARTUP <<EOF
from casatasks import casalog
casalog.setMemoryTotal(${mem_mb})
casalog.setMemoryFraction(99)
EOF
    echo "CASA cachedir $CASA_CACHEDIR: memory ${mem_mb} MB of ${alloc_mb} MB allocated, memfrac 99%"
}
make_casa_startup

run_casa () {
    xvfb-run -a $CASA_EXECUTABLE --cachedir $CASA_CACHEDIR --startupfile $CASA_STARTUP \
        --nologger --nogui --log2term -c "$@"
}

# Test imaging, sharded: the (field, config) tclean runs are independent, so run them in
# TEST_SHARDS parallel CASA processes. plan (once) -> N workers -> merge (once).
#  - Each worker has its OWN working directory (CASA writes casa-*.log and temp files to the cwd),
#    its OWN --cachedir and startup file, so the processes cannot clobber each other.
#  - Memory: each worker's CASA limit is 90% of the allocation / N; OpenMP threads = CPUS / N.
#  - Workers only write test_<field>_<config>.{psf,residual} (unique per job) into the shared
#    outdir; the merge step writes test_imaging_summary.txt. A crashed worker's jobs show up as
#    FAILED rows in the summary rather than stopping the job.
#  - Worker logs: <products>/test_imaging_shard_logs/shard_<i>.log
TEST_SHARDS=${TEST_SHARDS:-8}

run_test_shard () {      # <i> <n> <ms> <outdir> <shard_root>
    local i=$1 n=$2 ms=$3 out=$4 root=$5
    local wd=$root/shard_$i cache=$CASA_CACHEDIR/shard_$i
    mkdir -p $wd $cache
    cat > $cache/startup.py <<EOF
from casatasks import casalog
casalog.setMemoryTotal($(( ALLOC_MB * 9 / 10 / n )))
casalog.setMemoryFraction(99)
EOF
    local thr=$(( ${SLURM_CPUS_PER_TASK:-$n} / n ))
    [ $thr -ge 1 ] || thr=1
    cd $wd
    OMP_NUM_THREADS=$thr xvfb-run -a $CASA_EXECUTABLE --cachedir $cache --startupfile $cache/startup.py \
        --nologger --nogui --log2term -c $PIPE/run_test_imaging.py $ms $out shard $i/$n
}

run_test_imaging_sharded () {      # <ms> <outdir> <logdir>
    local ms=$1 out=$2 logdir=$3 n=$TEST_SHARDS
    if [ "$n" -le 1 ]; then
        run_casa $PIPE/run_test_imaging.py $ms $out; return
    fi
    local root=$WORKDIR/.test_shards_$(basename $ms .ms) pids=() i failed=0
    mkdir -p $out $logdir
    run_casa $PIPE/run_test_imaging.py $ms $out plan
    for (( i=0; i<n; i++ )); do
        ( run_test_shard $i $n $ms $out $root ) > $logdir/shard_$i.log 2>&1 &
        pids+=($!)
        sleep 3        # stagger: concurrent `xvfb-run -a` can race for the same display number
    done
    echo "$(date): $n test-imaging shards running (logs in $logdir)"
    for i in "${!pids[@]}"; do
        wait ${pids[$i]} || { echo "WARNING: test-imaging shard $i exited non-zero (see $logdir/shard_$i.log)"; failed=1; }
    done
    run_casa $PIPE/run_test_imaging.py $ms $out merge $n
    rm -rf $root
    [ $failed -eq 0 ] || echo "WARNING: some shards failed; affected jobs are marked FAILED in test_imaging_summary.txt"
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

    # Merge near-identical SPWs (sub-channel offsets) into one; no-op for 1 SPW.
    # Uses mstransform fftshift interpolation (fixed in CASA >= 6.7.5).
    run_casa $PIPE/combine_spws.py $WORKDIR/$MS_NAME

    run_casa $PIPE/find_weight_outliers.py $WORKDIR/$MS_NAME $PRODUCTS

    # Flagged % vs channel (channel number + frequency axes) after the hub-antenna flags
    run_casa $PIPE/plot_flag_vs_channel.py $WORKDIR/$MS_NAME $PRODUCTS \
        || echo "WARNING: plot_flag_vs_channel failed for $MS_NAME; continuing"

    if [[ "$MS_NAME" == *nocontsub* ]]; then
        echo "nocontsub MS: skipping uvcontsub"
    else
        run_casa $PIPE/run_uvcontsub.py $WORKDIR/$MS_NAME
    fi

    run_test_imaging_sharded $WORKDIR/$MS_NAME $PRODUCTS/test_images $PRODUCTS/test_imaging_shard_logs
done

# Summary report (markdown + PNG gallery per MS, overview.md/csv). On by default.
if [ "${MAKE_SUMMARY:-1}" = "1" ]; then
    echo "$(date): summary report"
    run_casa $PIPE/summarize_stage1_products.py $WORKDIR \
        || echo "WARNING: summary report failed; run summarize_stage1_products.py by hand"
    echo "Summary report: $WORKDIR/stage1_qa_reports/overview.md"
fi

echo "Job 1 complete. Inspect $WORKDIR/*_products, then run:"
echo "  sbatch $PIPE/job2_archive_and_cleanup.sh $TARGET $WORKDIR"
seff $SLURM_JOB_ID || true

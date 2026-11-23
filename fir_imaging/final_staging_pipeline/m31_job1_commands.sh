#!/bin/bash
# REFERENCE ONLY: the sbatch commands for the M31 C+D and A+B+C+D staging runs.
# Copy/paste them one at a time on the cluster. Running this file does nothing but print this note.
echo "m31_job1_commands.sh is a cheat sheet; copy the commands from it, do not run it."
exit 0

# Sizes are the sum of the per-config tars (A+B+C+D numbers are estimates until each concat finishes).
# "Peak scratch" ~ total + 2x the largest MS (uvcontsub +1x, combine_spws on hilores/OH up to +2x).
# All MSs stay in the workdir until job 2, so don't have every job below live at once: check
# `diskusage_report` and run a few at a time (job 2 frees the MS data).
# --time values marked (guess) are NOT measured; calibrate from `seff`/the log of the first run.
# Job 1 defaults: 24 h, 16 CPUs, 96 G, TEST_SHARDS=8. Filters are passed as env vars (inherited by sbatch);
# --job-name just makes the queue readable. Each job 1 needs its own job 2 (commands at the bottom).

PIPE=~/lglbs_hi_scripts/fir_imaging/final_staging_pipeline

# ---------------------------------------------------------------------------------------------------
# 0a. Make sure the cluster checkout has the NEW scripts (CONFIGS/LINES/DRY_RUN filters). An old copy
#     silently ignores these variables and runs every m31 tar (73 tars; this happened on 2026-11-22,
#     job 5319998), and DRY_RUN=1 on an old copy would start real work on the login node.
# ---------------------------------------------------------------------------------------------------
cd ~/lglbs_hi_scripts && git checkout main && git pull && git log --oneline -3 && cd -
grep -c DRY_RUN $PIPE/job1_prep_and_test.sh      # must print a number > 0; 0 means the old script

# ---------------------------------------------------------------------------------------------------
# 0b. Check the selection first (login node, no job submitted): tars, sizes, estimated scratch
# ---------------------------------------------------------------------------------------------------
DRY_RUN=1 CONFIGS=C+D bash $PIPE/job1_prep_and_test.sh m31
DRY_RUN=1 CONFIGS=A+B+C+D LINES=hilores bash $PIPE/job1_prep_and_test.sh m31

# Before running ABCD hilores: the tar was 314 G on 2026-11-22, but A+B+C+D inputs sum to 386 G.
# Make sure the concat has finished / the size stopped changing:
ls -l /home/ekoch/projects/rrg-eros-ab/ekoch/VLAXL/staged_measurement_sets/m31_A+B+C+D_*.ms.tar

# ---------------------------------------------------------------------------------------------------
# 1. C+D, all 8 lines in one job: 2.29 T, peak scratch ~3.1 T, 8 MSs
# ---------------------------------------------------------------------------------------------------
CONFIGS=C+D sbatch --job-name=m31_CD --time=72:00:00 $PIPE/job1_prep_and_test.sh m31      # (guess)

# If that is too long, split it (HI lines 1.18 T / the rest 1.11 T):
CONFIGS=C+D LINES=hi21cm,hi21cm_0p8kms sbatch --job-name=m31_CD_hi $PIPE/job1_prep_and_test.sh m31
CONFIGS=C+D LINES=hilores,himidres,oh1612,oh1665,oh1667,oh1720 sbatch --job-name=m31_CD_rest $PIPE/job1_prep_and_test.sh m31

# Smallest realistic first test (81 G): confirms sharded test imaging works on a real MS
CONFIGS=C+D LINES=hilores sbatch --job-name=m31_CD_hilores_test --time=12:00:00 $PIPE/job1_prep_and_test.sh m31   # (guess)

# ---------------------------------------------------------------------------------------------------
# 2. A+B+C+D, split by line (one concat tar per line; no nocontsub versions exist for the concats)
# ---------------------------------------------------------------------------------------------------
# hilores + himidres: ~1.1 T
CONFIGS=A+B+C+D LINES=hilores,himidres sbatch --job-name=m31_ABCD_hilo_himid $PIPE/job1_prep_and_test.sh m31

# OH lines: ~0.9-1.0 T each, peak scratch ~2.6-2.9 T
CONFIGS=A+B+C+D LINES=oh1612 sbatch --job-name=m31_ABCD_oh1612 $PIPE/job1_prep_and_test.sh m31
CONFIGS=A+B+C+D LINES=oh1665 sbatch --job-name=m31_ABCD_oh1665 $PIPE/job1_prep_and_test.sh m31
CONFIGS=A+B+C+D LINES=oh1667 sbatch --job-name=m31_ABCD_oh1667 $PIPE/job1_prep_and_test.sh m31
CONFIGS=A+B+C+D LINES=oh1720 sbatch --job-name=m31_ABCD_oh1720 $PIPE/job1_prep_and_test.sh m31

# HI, 0.8 km/s channels: ~1.6 T, peak scratch ~3.2 T
CONFIGS=A+B+C+D LINES=hi21cm_0p8kms sbatch --job-name=m31_ABCD_hi0p8 --time=48:00:00 $PIPE/job1_prep_and_test.sh m31   # (guess)

# HI full resolution: ~3.2 T, peak scratch ~6.4 T. Biggest job: check scratch quota first.
CONFIGS=A+B+C+D LINES=hi21cm sbatch --job-name=m31_ABCD_hi21cm --time=96:00:00 $PIPE/job1_prep_and_test.sh m31   # (guess)

# Optional overrides (any job): fewer/more test-imaging workers, more memory
#   TEST_SHARDS=4 ...                       (default 8; 1 = serial)
#   sbatch --mem=128G ...                   (if seff shows the 96 G is too tight)

# ---------------------------------------------------------------------------------------------------
# 3. Job 2: one per job 1, AFTER inspecting its workdir. Always pass the workdir explicitly.
# ---------------------------------------------------------------------------------------------------
# Workdirs are /home/ekoch/scratch/VLAXL_imaging/m31_<tag>_JOBID, e.g. m31_CD_JOBID, m31_ABCD_oh1612_JOBID.  (replace JOBID with the real job id)
# The pointer files hold them too:
cat /home/ekoch/scratch/VLAXL_imaging/last_workdir_m31_CD.txt
cat /home/ekoch/scratch/VLAXL_imaging/last_workdir_m31_ABCD_oh1612.txt
ls -d /home/ekoch/scratch/VLAXL_imaging/m31_*

sbatch --job-name=m31_CD_job2 --time=24:00:00 $PIPE/job2_archive_and_cleanup.sh m31 /home/ekoch/scratch/VLAXL_imaging/m31_CD_JOBID   # (guess; default 6 h)
sbatch --job-name=m31_ABCD_oh1612_job2 --time=24:00:00 $PIPE/job2_archive_and_cleanup.sh m31 /home/ekoch/scratch/VLAXL_imaging/m31_ABCD_oh1612_JOBID
# 3.2 T hi21cm: tar + verify + cp + cmp is several full passes over the data; 6 h is very likely too short
sbatch --job-name=m31_ABCD_hi21cm_job2 --time=72:00:00 $PIPE/job2_archive_and_cleanup.sh m31 /home/ekoch/scratch/VLAXL_imaging/m31_ABCD_hi21cm_JOBID   # (guess)

# Delete the workdir after job 2 succeeds (job 2 already removes the MSs and keeps the QA products):
PURGE=1 sbatch --export=ALL,PURGE=1 $PIPE/job2_archive_and_cleanup.sh m31 /home/ekoch/scratch/VLAXL_imaging/m31_CD_JOBID

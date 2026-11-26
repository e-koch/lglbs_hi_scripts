
## Staging from fir nearline

Long-term LGLBS storage of calibrated visibilities are on nearline at
`/home/ekoch/nearline/rrg-eros-ab/ekoch/calibrated`.

1. From the nearline folder (above), open the relevant target folder and run:

    ```
    for file in *.speclines.*.tar; do
        lfs hsm_restore $file;
    done
    ```

    Switch to `*.continuum.*.tar` for continuum data.

    Come back a few hr later and check in the same folder with `ls -lah` that all files have a reasonable size (>10 to 100s GB). See https://docs.alliancecan.ca/wiki/Using_nearline_storage for further information.

2. Transfer the target folder on nearline to `/home/ekoch/scratch/VLAXL_imaging/MeasurementSets`. These will tend to be several TB so globus is the preferred transfer method.

3. Run the staging portion of the pipeline:
 - Alter the line products to be created (if needed) in `run_lglbs_line_staging.py`
 - Submit the job passing the target name via the cmd line: `sbatch ~/lglbs_hi_scripts/fir_imaging/job_submit_lglbs_staging.sh wlm A+B+C+D` (NOTE: For wlm, ic1613, ic10, `ctr` will be appended by the script in the name to only select the central pointing with A and B config coverage).
 - Submit a separate job for the dwarfs that only have C+D hexagonal mosaic coverage (wlm, ic1613, ic10)


NOTE: `run_lglbs_line_staging.py` creates each line product separately (looping over multiple selections). This minimizes disc usage by cleaning up all intermediate products.

### M31

M31's LGLBS data volume exceeds the 20 TB scratch space. To mitigate this, run each config separataly,
then we will setup a manual script to concat them together.

1. Transfer M31 tracks per config to `scratch/VLAXL_imaging/MeasurementSets/m31`.

    Run the staging for the config:
    ```
    sbatch ~/lglbs_hi_scripts/fir_imaging/job_submit_lglbs_staging.sh m31 D
    ```

2. Split NGC185 and NGC205 out of the per-config tars. These two dwarfs were observed in the 14A-235 D-config tracks
   and must not be part of the M31 mosaic. The staging key (`ms_file_key.txt`) now selects `M31*` only, so D tars staged
   after that fix are already clean; the D tars staged earlier still contain them (and so do any concats built from them).
   Run this once the D config is staged, before step 4:

    ```
    sbatch ~/lglbs_hi_scripts/fir_imaging/job_split_m31_ngc_fields.sh
    ```

    - Processes `m31_D_<line>[_nocontsub].ms.tar` for every line (the concatenated `m31_C+D_*` / `m31_A+B+C+D_*` tars are never touched). For each tar it extracts the MS, `split`s the `M31*` fields and each NGC field into their own MS (all columns, flags kept, no rows dropped), checks that the row counts add up and no MAIN column was lost, then writes `ngc185_D_<line>[_nocontsub].ms.tar` and `ngc205_D_<line>[_nocontsub].ms.tar` into `staged_measurement_sets` and replaces the `m31_D_*` tar with the `M31*`-only one (last, so an interrupted job can be repeated).
    - Safe to repeat: a tar with no NGC fields with data is skipped.
    - To check A, B and C for NGC fields first (only reads each tar's FIELD table, writes nothing; the FIELD table may list fields without data, which the real run would skip):

    ```
    M31_SPLIT_DRY_RUN=1 M31_SPLIT_CONFIGS=A,B,C,D sbatch --time=02:00:00 ~/lglbs_hi_scripts/fir_imaging/job_split_m31_ngc_fields.sh
    ```

    - Env vars: `M31_SPLIT_CONFIGS` (default `D`), `M31_SPLIT_LINES` (default all), `M31_SPLIT_DRY_RUN`, `M31_SPLIT_STAGED_DIR` (test copy of the staging dir). Work happens in `scratch/VLAXL_imaging/m31_ngc_split_<jobid>` (MSs are removed by the script; the dir with the CASA logs is left behind). Peak scratch is about 2.2x the largest D MS. Same CASA (6.7.6) setup as the concat job.
    - The new `ngc185_*`/`ngc205_*` tars are named like any other target's, so `job1_prep_and_test.sh ngc185` (or `ngc205`) will pick them up; the `m31` runs will not.
    - After this, rerun the concat jobs (step 4) and job 1 for the concatenated tars. The existing `m31_C+D_*` / `m31_A+B+C+D_*` tars and any job 1 workdir built from them (e.g. `m31_CD_6016625`) contain the NGC fields: do not run job 2 on them.

3. Iterate through all configs.

4. Concatenate individual configs together. Submit one job per combination (they can run in parallel):

    ```
    sbatch ~/lglbs_hi_scripts/fir_imaging/job_submit_m31_concat.sh ABCD
    sbatch ~/lglbs_hi_scripts/fir_imaging/job_submit_m31_concat.sh CD
    ```

    - `ABCD` creates `m31_A+B+C+D_<line>.ms.tar`; `CD` creates `m31_C+D_<line>.ms.tar`.
    - For each line product, the per-config tars are copied from `staged_measurement_sets`, extracted, concatenated, re-tarred, and the combined tar is moved back to `staged_measurement_sets`.
    - Each job works in its own scratch dir, `scratch/VLAXL_imaging/m31_concat_<combo>_<jobid>`. The MSs are cleaned up by the script, but the dir (with CASA logs) is left behind; remove it by hand when finished.
    - The combo is passed to `run_m31_concat_configs.py` through the `M31_CONCAT_COMBO` and `M31_CONCAT_LOCAL_DIR` environment variables, set by the submit script.
    - Uses the same CASA (6.7.6) and startup setup as `final_staging_pipeline/job1_prep_and_test.sh`.

    TODO: run a final statwt.


## Example submission for dwarf target

sbatch  ~/lglbs_hi_scripts/fir_imaging/job_submit_lglbs_staging.sh ic1613 A+B+C+D
sbatch --time=12:00:00 --mem=32G ~/lglbs_hi_scripts/fir_imaging/job_submit_lglbs_staging.sh ic1613 C+D


# Imaging on fir


# Test imaging on fir

Submit the test imaging job, providing the target, config, and line name as defined in the keys.


sbatch  ~/lglbs_hi_scripts/fir_imaging/job_submit_lglbs_HI_imaging_param_test.sh ngc6822 A+B+C+D hilores


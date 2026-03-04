
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

2. Iterate through all configs.

3. Concatenate individual configs together.

    TODO: setup script to do this an run a final statwt.


## Example submission for dwarf target

sbatch  ~/lglbs_hi_scripts/fir_imaging/job_submit_lglbs_staging.sh ic1613 A+B+C+D
sbatch --time=12:00:00 --mem=32G ~/lglbs_hi_scripts/fir_imaging/job_submit_lglbs_staging.sh ic1613 C+D


# Imaging on fir


# Test imaging on fir

Submit the test imaging job, providing the target, config, and line name as defined in the keys.


sbatch  ~/lglbs_hi_scripts/fir_imaging/job_submit_lglbs_HI_imaging_param_test.sh ngc6822 A+B+C+D hilores


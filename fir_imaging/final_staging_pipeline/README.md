# Final staging pipeline

Post-staging QA for the tars in `staged_measurement_sets`. Two jobs so job-1
outputs can be checked before the archive is overwritten.

    sbatch job1_prep_and_test.sh <target>
    # inspect /home/ekoch/scratch/VLAXL_imaging/<target>_<jobid>/*_products/
    sbatch job2_archive_and_cleanup.sh <target> [workdir]   # PURGE=1 to delete workdir

Job 1, per staged MS tar: extract to the work dir -> `find_weight_outliers.py`
(flag hub antennas, `*_flag_summary.txt` before/after) -> `run_uvcontsub.py`
(order 1, first/last 5% of channels per SPW; skipped for `nocontsub`) ->
`run_test_imaging.py` (natural-weighted single-channel PSF + residual per
field per VLA config; configs by absolute p95 baseline).
Job 2: tar, verify, replace tar in project space, remove local MS.

`summary_amp_weight_vs_uvdist.py` is kept here as a standalone diagnostic.
Set `STAGED_DIR` to point both jobs at a test copy of the archive.

`summarize_stage1_products.py [DATA_DIR] [-o OUTDIR]` builds a per-run markdown summary (flagging
before/after, hub-antenna flag commands, test-imaging table, attention flags) and a PNG gallery of
the PSF/residual test images for every `*_products` dir, plus `overview.md`/`overview.csv`. It reads
CASA images with casatools if available, otherwise casa-formats-io.

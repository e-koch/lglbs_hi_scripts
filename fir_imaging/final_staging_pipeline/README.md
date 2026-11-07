# Final staging pipeline

Post-staging QA for the tars in `staged_measurement_sets`. Two jobs so job-1
outputs can be checked before the archive is overwritten.

    sbatch job1_prep_and_test.sh <target>
    # inspect /home/ekoch/scratch/VLAXL_imaging/<target>_<jobid>/*_products/
    sbatch job2_archive_and_cleanup.sh <target> [workdir]   # PURGE=1 to delete workdir

Job 1 runs CASA 6.7.6 (`casa-6.7.6-14-py3.12.el9`). Per staged MS tar: extract to the work dir ->
`combine_spws.py` (if >1 SPW is in use, e.g. OH lines with a few B0D0 scans ~0.25 chan from the A0C0 SPW,
or M31 hilores; mstransform regrid onto the busiest SPW's grid with `fftshift` interpolation, in place;
stops without changes if SPWs differ in nchan/width or are offset by >=1 chan) -> `find_weight_outliers.py`
(flag hub antennas, `*_flag_summary.txt` before/after) -> `plot_flag_vs_channel.py`
(flagged % vs channel number + frequency, PNG/CSV per SPW) -> `run_uvcontsub.py`
(order 1, first/last 5% of channels per SPW; skipped for `nocontsub`) ->
`run_test_imaging.py` (natural-weighted single-channel PSF + residual per
field per VLA config; configs by absolute p95 baseline).
Job 2: tar, verify, replace tar in project space, remove local MS.

`summary_amp_weight_vs_uvdist.py` is kept here as a standalone diagnostic.
Set `STAGED_DIR` to point both jobs at a test copy of the archive.

Job 1 runs `summarize_stage1_products.py` on its work dir at the end (`MAKE_SUMMARY=0` to skip).
`summarize_stage1_products.py [DATA_DIR] [-o OUTDIR]` builds a per-run markdown summary (flagging
before/after, hub-antenna flag commands, test-imaging table, attention flags) and a PNG gallery of
the PSF/residual test images for every `*_products` dir, plus `overview.md`/`overview.csv`. It reads
CASA images with casatools if available, otherwise casa-formats-io.

# Final staging pipeline

Post-staging QA for the tars in `staged_measurement_sets`. Two jobs so job-1
outputs can be checked before the archive is overwritten.

    sbatch job1_prep_and_test.sh <target>
    # inspect /home/ekoch/scratch/VLAXL_imaging/<target>_<jobid>/*_products/
    sbatch job2_archive_and_cleanup.sh <target> [workdir]   # PURGE=1 to delete workdir

Job 1 runs CASA 6.7.6 (`casa-6.7.6-14-py3.12.el9`). Per staged MS tar: extract to the work dir ->
`combine_spws.py` (if >1 SPW is in use, e.g. OH lines with a few B0D0 scans ~0.25 chan from the A0C0 SPW,
or M31 hilores; mstransform regrid onto the busiest SPW's grid with `fftshift` interpolation, in place;
far-offset SPWs (>=1 chan, or different nchan/width) that are recorded in the same scans as the busiest SPW
are dropped with `split` instead of merged, e.g. the neighbouring OH1667 SPW that the OH1665 window picks up in
M33 17B-162; stops without changes if an offset SPW is not simultaneous with the busiest SPW) -> `find_weight_outliers.py`
(flag hub antennas, `*_flag_summary.txt` before/after) -> `plot_flag_vs_channel.py`
(flagged % vs channel number + frequency, PNG/CSV per SPW) -> `run_uvcontsub.py`
(order 1, first/last 5% of channels per SPW; skipped for `nocontsub`) ->
`run_test_imaging.py` (natural-weighted single-channel PSF + residual per
field per VLA config; configs by absolute p95 baseline).
Job 2: tar, verify, replace tar in project space, remove local MS, then copy the QA products to
`staged_measurement_sets_summary/<target>/` (`*_products/` minus `test_images/`, plus `stage1_qa_reports/`;
override the root with `SUMMARY_DIR`).

Large targets (M31): job 1 runs every `<target>_*.ms.tar` by default. Select a subset with env vars
(`CONFIGS=C+D`, `LINES=hilores,himidres`, `CONTSUB=yes|no`; `DRY_RUN=1 bash job1_prep_and_test.sh m31`
lists the selection, sizes and estimated scratch need). Each submission gets its own workdir
(`<target>_<tag>_<jobid>`) and needs its own job 2, passing the workdir explicitly. All MSs stay in
the workdir until job 2, so scratch use is ~ the sum of the selected tars plus 1-2x the current MS.
e.g. `CONFIGS=A+B+C+D LINES=oh1612 sbatch job1_prep_and_test.sh m31`.

Test imaging is sharded over `TEST_SHARDS` (default 8) CASA processes, each with its own working
dir, CASA cache dir and 1/N of the memory budget: `run_test_imaging.py <ms> <outdir> plan`, then
`shard I/N` x N in parallel, then `merge N` (no mode = serial, as in `test_imaging_single_ms.sh`).
Worker logs go to `<name>_products/test_imaging_shard_logs/`.

`summary_amp_weight_vs_uvdist.py` is kept here as a standalone diagnostic.
Set `STAGED_DIR` to point both jobs at a test copy of the archive.

Job 1 runs `summarize_stage1_products.py` on its work dir at the end (`MAKE_SUMMARY=0` to skip).
`summarize_stage1_products.py [DATA_DIR] [-o OUTDIR]` builds a per-run markdown summary (flagging
before/after, hub-antenna flag commands, test-imaging table, attention flags) and a PNG gallery of
the PSF/residual test images for every `*_products` dir, plus `overview.md`/`overview.csv`. It reads
CASA images with casatools if available, otherwise casa-formats-io.

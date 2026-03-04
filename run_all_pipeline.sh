
# Need to be on uvcontsub_exclude_freqrange for now
git clone https://github.com/e-koch/phangs_imaging_scripts.git
git checkout origin/uvcontsub_exclude_freqrange
git switch -c uvcontsub_exclude_freqrange


/opt/casa-6.2.1-7-pipeline-2021.2.0.128/bin/casa --nogui --log2term -c ~/lglbs_hi_scripts/run_lglbs_staging_dwarfctr.py ic10ctr
/opt/casa-6.2.1-7-pipeline-2021.2.0.128/bin/casa --nogui --log2term -c ~/lglbs_hi_scripts/run_lglbs_staging_dwarfctr.py ic1613ctr
/opt/casa-6.2.1-7-pipeline-2021.2.0.128/bin/casa --nogui --log2term -c ~/lglbs_hi_scripts/run_lglbs_staging_dwarfctr.py wlmctr
# ngc6822 has no centre only pointing


# Derived products for feathered C+D full dwarf mosaics
conda activate pyuvdata_ewk
python ~/lglbs_hi_scripts/run_lglbs_derived.py ic10 &
python ~/lglbs_hi_scripts/run_lglbs_derived.py ic1613 &
python ~/lglbs_hi_scripts/run_lglbs_derived.py ngc6822 &
python ~/lglbs_hi_scripts/run_lglbs_derived.py wlm &


# Make tar files for the C+D full dwarf mosaics
# And upload to the archive
cd /reduction10/erickoch/LGLBS/line_imaging/derived
echo "Making tar files"
tar -cf ic10_c+d+tp_hi.tar ic10
tar -cf ic1613_c+d+tp_hi.tar ic1613
tar -cf ngc6822_c+d+tp_hi.tar ngc6822
tar -cf wlm_c+d+tp_hi.tar wlm

~/rclone-v1.65.2-linux-amd64/rclone copy ic10_c+d+tp_hi.tar lglbs-gdrive:"LineImaging/HI_C+D_feathered/" --progress
~/rclone-v1.65.2-linux-amd64/rclone copy ic1613_c+d+tp_hi.tar lglbs-gdrive:"LineImaging/HI_C+D_feathered/" --progress
~/rclone-v1.65.2-linux-amd64/rclone copy ngc6822_c+d+tp_hi.tar lglbs-gdrive:"LineImaging/HI_C+D_feathered/" --progress
~/rclone-v1.65.2-linux-amd64/rclone copy wlm_c+d+tp_hi.tar lglbs-gdrive:"LineImaging/HI_C+D_feathered/" --progress

# And upload the directories
~/rclone-v1.65.2-linux-amd64/rclone copy ic10 lglbs-gdrive:"LineImaging/HI_C+D_feathered/ic10/" --progress
~/rclone-v1.65.2-linux-amd64/rclone copy ic1613 lglbs-gdrive:"LineImaging/HI_C+D_feathered/ic1613/" --progress
~/rclone-v1.65.2-linux-amd64/rclone copy ngc6822 lglbs-gdrive:"LineImaging/HI_C+D_feathered/ngc6822/" --progress
~/rclone-v1.65.2-linux-amd64/rclone copy wlm lglbs-gdrive:"LineImaging/HI_C+D_feathered/wlm/" --progress

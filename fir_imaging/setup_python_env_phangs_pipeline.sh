
cd

# python 3.13 install
module load StdEnv/2023

module load python/3.13

# ipython/numpy/scipy/mpl
module load scipy-stack


virtualenv --no-download venv_phangs_pipeline
source $HOME/venv_phangs_pipeline/bin/activate

# Check casa packages are available
avail_wheels 'casa*'

pip install --no-index --upgrade pip

# (09/2026) Until we push a new release of casa-formats-io, it needs to be compiled.
# git clone https://github.com/radio-astro-tools/casa-formats-io.git
cd casa-formats-io
# git pull origin main
pip install .
cd ../

# Install from phangs-pipeline
cd phangs_imaging_scripts
# git pull origin master
pip install -e .[casa]
cd ../


# Check import works
python -c 'import casatasks'
python -c 'import spectral_cube'
python -c 'import phangsPipeline'
# Works!

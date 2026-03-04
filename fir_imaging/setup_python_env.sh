
cd

# python 3.11 install
module load StdEnv/2023

# ipython/numpy/scipy/mpl
module load scipy-stack


virtualenv --no-download phangs
source $HOME/phangs/bin/activate

pip install --no-index --upgrade pip

# astropy=7.2 needs numpy 2.4
pip install 'astropy<7.2'

pip install spectral-cube

# Check import works
python -c 'import spectral_cube'

# Issue is the build of casa-formats-io is old
# compiled with numpy-1.X.
git clone https://github.com/radio-astro-tools/casa-formats-io.git
cd casa-formats-io
pip install .

# Works!

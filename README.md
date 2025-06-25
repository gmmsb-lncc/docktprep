# docktprep
Prepare and fix [Protein Data Bank (PDB)](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/beginner%E2%80%99s-guide-to-pdbx-mmcif) files for molecular modeling and docking.

## Installation
Clone this repository, create a virtual environment and install the required dependencies:
```bash
# clone this repository...
cd docktprep
python3.12 -m venv env  # use python >=3.10 
source env/bin/activate
python -m pip install -r requirements.txt
```

Run tests:
```bash
python -m pytest -vs tests/ --log-cli-level=INFO 
```

Run the app:
```bash
python -m docktprep.main --help
```


### MODELLER installation
Some functionalities of this package depend on [MODELLER](https://salilab.org/modeller/), which is not included in this repository; you'll need to install it separately. Follow the [installation instructions](https://salilab.org/modeller/download_installation.html) on the official MODELLER website.

> [!TIP]
> Below, instructions for Debian/Ubuntu systems are provided:
> ```bash
> wget https://salilab.org/modeller/10.7/modeller_10.7-1_amd64.deb
> sudo env KEY_MODELLER=XXXX dpkg -i modeller_10.7-1_amd64.deb  # substitute 'XXXX' for the actual key
> python3.12 -c "import modeller; print(modeller.__version__)"  # test the installation
> ```

> [!IMPORTANT]
> Once MODELLER is installed, you must create symbolic links to its files in the virtual environment’s site-packages directory to make the library available inside the venv.

Find the location of the MODELLER installation:
```bash
deactivate  # `deactivate` the virtual environment first, if it is active!
MODELLER_PATH=$(python3.12 -c "import modeller, os; print(os.path.dirname(modeller.__file__))")
MODELLER_PARENT=$(dirname -- "${MODELLER_PATH%/}")
```

Reactivate the virtual environment and create symbolic links:
```bash
source env/bin/activate
VENV_SITEPACKAGES=$(python -c "import site; print(site.getsitepackages()[0])")
ln -s "$MODELLER_PATH" "$VENV_SITEPACKAGES/modeller"
ln -s "$MODELLER_PARENT/_modeller.so" "$VENV_SITEPACKAGES/_modeller.so"
```

Test the installation inside the virtual environment:
```bash
python -c "import modeller; print(modeller.__version__)"
```








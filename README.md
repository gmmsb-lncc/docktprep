# docktprep
> ⚠️ Under development.

Create DockThor input files from PDB, mmCIF and MOL2 formats. Prepare and fix common issues in the PDB files.

## Installation
Clone this repository:
```
git clone https://github.com/gmmsb-lncc/docktprep.git
```

Create a virtual environment and install the required dependencies:
```
cd docktprep
python3 -m venv env
source env/bin/activate
python -m pip install -r requirements.txt
```

Run tests:
```
python -m pytest -vs tests/ --log-cli-level=INFO 
```

Run the app:
```
python -m docktprep.main --help
```



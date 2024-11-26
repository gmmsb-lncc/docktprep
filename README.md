# docktprep

> This is a work in progress and the code is not yet stable.
  
 Create DockThor input files from PDB, mmCIF and MOL2 formats. Prepare and fix common issues in the PDB files.

 Accepted input formats:
 - Receptor: PDB, PDBx/mmCIF ⏳.
 - Ligand: PDB ⏳, PDBx/mmCIF ⏳, MOL2 ⏳.


Features:
- Handles disordered atom positions.¹
- Handles amino acid residue insertions.¹
- Handles non-standard amino acids nomenclatures (e.g. HISD, HSD). ⏳
- Allows for structure model selection.¹
- Supports DNA/RNA receptors.¹
- Replaces non-standard amino acids/nucleotides with standard ones.²
- Adds missing heavy atoms to the structure.²
- Adds missing residues to the structure.² ⏳
- Adds missing hydrogens to the structure.²

Features originally implemented in ¹[Biopython](https://biopython.org/) and ²[PDBFixer](https://github.com/openmm/pdbfixer).

## Installation and usage as a command-line tool

To install the application as a command-line tool, you'll need [conda](https://docs.anaconda.com/miniconda/install/). 

First, create a new conda environment from the `environment.yml` file:
```bash
conda env create --file environment.yml --prefix ./venv/
```

You will always need to activate the environment before running the application:
```bash
conda activate ./venv/
```

That's it. For a complete list of options, run:
```bash
docktprep --help
```

## Development setup
Follow the instructions below to set up the development environment.

### Installation
Clone the repository and install the dependencies using [conda](https://docs.anaconda.com/miniconda/install/).

```bash
# create a new conda environment from the environment.yml file
conda env create --file environment.yml --prefix ./venv/
conda activate ./venv/

```

### Run tests

```bash
python -m pytest -vs tests/ --log-cli-level=INFO
```

### Execute main script

```bash
python -m docktprep.main --help
```

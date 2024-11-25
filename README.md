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

## Installation

```bash
python -m pip install docktprep
```

## Usage

For a complete list of options, run:
```bash
docktprep --help
```

## Development setup

### Installation
Clone the repository and install the dependencies using [conda](https://docs.anaconda.com/miniconda/install/).

```bash
# create a new conda environment
conda create --prefix ./venv/ python=3.12
conda activate ./venv/

# install the dependencies from .yml file
conda env update --file environment.yml
```

### Run tests

```bash
python -m pytest -vs tests/ --log-cli-level=INFO
```

### Execute main script

```bash
python -m docktprep.main --help
```

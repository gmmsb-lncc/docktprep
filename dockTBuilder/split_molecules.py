import os
import sys
import pandas as pd
import biopandas as bdp
from biopandas.pdb import PandasPdb
from rdkit import Chem
from rdkit.Chem import AllChem

def pdb_exist(pathName):
    return os.path.isfile(pathName)

def split_file(pathName):
    ppdb = PandasPdb()
    ppdb.read_pdb(pathName)

    atoms = ppdb.df['ATOM']
    hetatms = ppdb.df['HETATM']

    hetatms_split = hetatms.groupby("residue_name")
                        
    return atoms, hetatms_split

def add_hydrogens(hetatms_df):

    mol_block = PandasPdb().to_pdb(hetatms_df)
    mol = Chem.MolFromPDBBlock(mol_block)

    mol_with_hydrogens = AllChem.AddHs(mol, addCoords=True, addResidueInfo=True)

    # Converter o molécula com hidrogênios de volta para um DataFrame
    hetatms_with_df = PandasPdb().pdb_to_df(mol_with_hydrogens)

    return hetatms_with_hydrogens_df

def save_mol_as_pdb(mol, output_file):
    pdb_block = Chem.MolToPDBBlock(mol)

    with open(output_file, 'w') as file:
        file.write(pdb_block)

def add_hydrogens(hetatms_df):
    ppdb = PandasPdb()
    ppdb.df['HETATM'] = hetatms_df
    tmp_file: str = "tmp.pdb"

    ppdb.to_pdb(path=tmp_file, records=None, gz=False, append_newline=True)

    with open(tmp_file, "r") as f:
        file_content = f.read()
        mol = Chem.MolFromPDBBlock(file_content)

    mol_with_hydrogens = Chem.rdmolops.AddHs(mol, addCoords=True, addResidueInfo=True)
    os.remove(tmp_file)

    return mol_with_hydrogens

def main():
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <path_to_pdb_file>")
        sys.exit(1)

    moleculePath = sys.argv[1]

    if pdb_exist(moleculePath):
        print(f"The {moleculePath} exist!\n")
    else:
        print(f"The {moleculePath} doesn't exist! Try it again inserting the correct file path.\n")
        sys.exit(1)

    atoms_df, hetatms_df = split_file(moleculePath)

    for residue_name, group_hetatm in hetatms_df:
        if group_hetatm.empty:
            continue
        
        mol_with_hydrogens = add_hydrogens(group_hetatm)
        save_mol_as_pdb(mol_with_hydrogens, f"{residue_name}_hs.pdb")
    


if __name__ == "__main__":
    main()
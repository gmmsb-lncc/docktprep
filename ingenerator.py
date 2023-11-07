import sys
from prody import *
from rdkit import Chem
from rdkit.Chem import AllChem
from io import StringIO
import pypdb


macromolecules = {
    'amino_acids': {
        'ALA': 'Alanine',
        'ARG': 'Arginine',
        'ASN': 'Asparagine',
        'ASP': 'Aspartic Acid',
        'CYS': 'Cysteine',
        'GLN': 'Glutamine',
        'GLU': 'Glutamic Acid',
        'GLY': 'Glycine',
        'HIS': 'Histidine',
        'ILE': 'Isoleucine',
        'LEU': 'Leucine',
        'LYS': 'Lysine',
        'MET': 'Methionine',
        'PHE': 'Phenylalanine',
        'PRO': 'Proline',
        'SER': 'Serine',
        'THR': 'Threonine',
        'TRP': 'Tryptophan',
        'TYR': 'Tyrosine',
        'VAL': 'Valine'
    },
    'nucleic_acids': {
        'A': 'Adenine',
        'C': 'Cytosine',
        'G': 'Guanine',
        'T': 'Thymine',
        'U': 'Uracil'
    },
    'carbohydrates': {
        'GLC': 'Glucose',
        'MAN': 'Mannose',
        'GAL': 'Galactose',
        'FUC': 'Fucose',
        'NAG': 'N-Acetylglucosamine',
        'NDG': 'N-Diacetylglucosamine',
        'BMA': 'Beta-D-Mannose',
        'SIA': 'Sialic Acid'
    }
}

def get_macromolecule_type(resname):
    resname_upper = resname.upper()
    if resname_upper in macromolecules['amino_acids']:
        return 'Aminoácido'
    elif resname_upper in macromolecules['nucleic_acids']:
        return 'Ácido Nucleico'
    elif resname_upper in macromolecules['carbohydrates']:
        return 'Carboidrato'
    else:
        return 'Desconhecido'
        
def listMolecules(moleculeName):

    parser = PDB.PDBParser(QUIET=True)

    try:
        structure = parser.get_structure("structure", moleculeName)

        receptors = {}  
        ligands = set()


        for atom in structure.get_atoms():
            if atom.element != "H":
                resname = atom.parent.resname
                if resname.strip():
                    if is_macromolecule(resname):
                        macromolecule_type = get_macromolecule_type(resname)
                        if resname not in receptors:
                            receptors[resname] = macromolecule_type
                    else:
                        ligands.add(resname)

        print("Moléculas do Receptor:")
        for resname, macromolecule_type in receptors.items():
            print(f"{resname} - {macromolecule_type}")

        print("\nLigantes:")
        for ligand in ligands:
            print(ligand)

    except FileNotFoundError:
        print(f'O arquivo PDB "{moleculeName}" não foi encontrado.')


def add_hydrogens(molecule):
    Chem.AddHs(molecule)
    return molecule


if __name__ == "__main__":
    moleculeName = sys.argv[1]        
    #listMolecules(moleculeName)

    parser = PDB.PDBParser(QUIET=True)

    try:
        structure = parser.get_structure("structure", moleculeName)
    
        # Itera sobre os modelos, cadeias e átomos e adiciona hidrogênios usando RDKit
        for model in structure:
            for chain in model:
                for residue in chain:
                    residue_name = residue.resname
                    residue_id = residue.id[1]
                    pdb_id = f"{residue_name} {residue_id}"

                    rdkit_molecule = Chem.MolFromPDBBlock(pdb_id, sanitize=False)

                    if rdkit_molecule is not None:
                        print("entra aqui?")
                        # Adiciona hidrogênios à molécula
                        rdkit_molecule_with_h = add_hydrogens(rdkit_molecule)
                        
                        # Atualiza as coordenadas dos átomos no arquivo PDB com hidrogênios
                        for atom, new_atom in zip(residue, rdkit_molecule_with_h.GetAtoms()):
                            print("entra aqui?")
                            atom.set_coord((new_atom.GetIdx(),) + tuple(new_atom.GetPos()))

            io = PDB.PDBIO()
            io.set_structure(structure)
            io.save("nova.pdb")

    except FileNotFoundError:
        print(f'O arquivo PDB "{moleculeName}" não foi encontrado.')            

                
    #             if rdkit_molecule is not None:
    #                 # Adiciona hidrogênios à molécula
    #                 rdkit_molecule_with_h = add_hydrogens(rdkit_molecule)
                    
    #                 # Atualiza as coordenadas dos átomos no arquivo PDB com hidrogênios
    #                 for atom, new_atom in zip(residue, rdkit_molecule_with_h.GetAtoms()):
    #                     atom.set_coord((new_atom.GetIdx(),) + tuple(new_atom.GetPos()))
                        
    # # Salva a estrutura PDB com hidrogênios
    # io = PDB.PDBIO()
    # io.set_structure(structure)
    # io.save(moleculeName_H.pdb)
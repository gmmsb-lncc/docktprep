import os
import requests
import sys
from Bio.PDB import PDBParser
from io import StringIO

def pdb_exist(pathName):
    return os.path.isfile(pathName)

def download_pdb(pdbCode):
    url = f"https://files.rcsb.org/download/{pdbCode}.pdb"
    try:
        response = requests.get(url)
        response.raise_for_status()
        with open(f"{pdbCode}.pdb", "wb") as pdb_file:
            pdb_file.write(response.content)
        print(f"The {pdbCode}.pdb was successfully downloaded!")
    except requests.exceptions.HTTPError as errh:
        print(f"HTTP error: {errh}")
        sys.exit()
    except requests.exceptions.ConnectionError as errc:
        print(f"Connection Error: {errc}")
        sys.exit()
    except requests.exceptions.Timeout as errt:
        print(f"Timeout Request: {errt}")
        sys.exit()
    except requests.exceptions.RequestException as err:
        print(f"Request error: {err}")
        sys.exit()

def split_file(pdbFile):
    atoms = []
    hetatms = {}
    
    with open(pdbFile, 'r') as pdbFile:
        for line in pdbFile:
            if line.startswith('ATOM'):
                atoms.append(line)
            elif line.startswith('HETATM'):
                hetatm_id = line[17:20].strip()  # Obtém o ID do HETATM
                if hetatm_id not in hetatms:
                    hetatms[hetatm_id] = []
                hetatms[hetatm_id].append(line)
                        
    return atoms, hetatms

def create_mol(pdb_string):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("mol", StringIO(pdb_string))
    return structure

moleculePath = sys.argv[1]

if pdb_exist(moleculePath):
    print(f"The {moleculePath} exists!\n")
else:
    answer = input(f"The {moleculePath} doesn't exist! Do you want to download {moleculePath} from RCSB PDB? (S/N)\n").strip().lower()
    if answer == 's':
        pdbCode = input(f"Enter the PDB code: ").strip().lower()
        download_pdb(pdbCode)
        moleculePath = f"{pdbCode}.pdb"
    else:
        print("Bye")
        sys.exit()

atoms, hetatms_dict = split_file(moleculePath)

# Parse dos átomos em um objeto PDB
atomsToPDB = "\n".join(atoms)
macromolecule = create_mol(atomsToPDB)

# Parse dos heteroátomos em objetos PDB com base no ID
smallMolecules = {}
for hetatmID, hetatmLines in hetatms_dict.items():
    hetatmToPDB = "\n".join(hetatmLines)
    smallMolecules[hetatmID] = create_mol(hetatmToPDB)

    
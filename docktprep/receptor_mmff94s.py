import pickle
import os
import io
import warnings
import logging

from Bio.PDB import PDBExceptions

from docktprep.receptor_parser import Receptor

def load_mmff94s_dict():
    """
    Load the MMFF94S dictionary from a pickle file.
    
    Returns:
        dict: A dictionary containing MMFF94S parameters.
    """
    with open("docktprep/mmff94s_dict.pkl", "rb") as f:
        mmff94s_dict = pickle.load(f)
    return mmff94s_dict

def write_topology(receptor: Receptor, file: str):
    """
    Write the topology file for the receptor using MMFF94S parameters.
    
    Args:
        receptor (Receptor): The receptor object containing the structure.
    """
    mmff94s_dict = load_mmff94s_dict()

    file_id = os.path.splitext(os.path.basename(receptor.file))[0] + "_topology"
    receptor.current_file_stream.seek(0)

    parser = receptor.get_biopython_parser()
    with warnings.catch_warnings(record=True) as warns:
        warnings.simplefilter("always")
        if receptor.file_ext == ".cif":
            structure = parser.get_structure(file_id, receptor.current_file_stream)
            mmcif_dict = parser._mmcif_dict.copy()
        elif receptor.file_ext == ".pdb":
            # Parse PDB as mmCIF by converting to mmCIF format in-memory
            pdb_content = receptor.current_file_stream.read()
            receptor.current_file_stream.seek(0)

            # Use Bio.PDB to parse PDB, then write to mmCIF in-memory
            pdb_parser = parser  # assuming parser is PDBParser
            structure = pdb_parser.get_structure(file_id, io.StringIO(pdb_content))

            # Write structure to mmCIF in-memory
            cif_stream = io.StringIO()
            mmcif_io = receptor.get_biopython_file_io(output_fmt="cif")
            mmcif_io.set_structure(structure)
            mmcif_io.save(cif_stream)
            cif_stream.seek(0)

            # Now parse mmCIF from in-memory stream
            mmcif_parser = receptor.get_biopython_parser(fmt="cif")
            structure = mmcif_parser.get_structure(file_id, cif_stream)
            mmcif_dict = mmcif_parser._mmcif_dict.copy()

    for warn in warns:
        if warn.category == PDBExceptions.PDBConstructionWarning:
            logging.warning(f"{receptor.file}: {str(warn.message).split("\n")[0]}")

    MMFF_atom_type = []
    MMFF_partial_charge = []
    for residue in structure.get_residues():
        res_name = residue.get_resname()
        if res_name in mmff94s_dict:
            for atom in residue.get_atoms():
                atom_name = atom.get_name()
                if atom_name in mmff94s_dict[res_name]:
                    atom_type, partial_charge = mmff94s_dict[res_name][atom_name]
                    MMFF_atom_type.append(atom_type)
                    MMFF_partial_charge.append(partial_charge)
    
    mmcif_dict['_atom_site.MMFF_atom_type'] = MMFF_atom_type
    mmcif_dict['_atom_site.MMFF_partial_charge'] = MMFF_partial_charge

    mmcif_io = receptor.get_biopython_file_io(output_fmt="cif")
    mmcif_io.set_dict(mmcif_dict)
    with open(file.split(".")[0] + "_topology.cif", "w") as f:
        mmcif_io.save(f)
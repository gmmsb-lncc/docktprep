import io
import logging
import os
import warnings
from typing import Protocol

from Bio.PDB import PDBExceptions
from modeller import *
from modeller.scripts import complete_pdb

from . import nonstd_residues
from .receptor_parser import FileFormatHandler, Receptor
from .stdout_manager import capture_output, suppress_output

__all__ = [
    "ModellerOperation",
    "CompletePDBOperation",
    "AddMissingAtomsOperation",
    "ReplaceNonStdResiduesOperation",
]


class ModellerOperation(Protocol):
    def run_modeller(self, receptor: Receptor, **kwargs) -> None: ...


class CompletePDBOperation(ModellerOperation):
    """Complete the PDB file by adding missing atoms and residues."""

    def run_modeller(self, receptor: Receptor, transfer_res_num: bool = False) -> None:
        with suppress_output():
            env = Environ()
            env.libs.topology.read(file="$(LIB)/top_heav.lib")
            env.libs.topology.read(file="$(LIB)/top_allh.lib")
            env.libs.parameters.read(file="$(LIB)/par.lib")
            env.io.hetatm = True
            env.io.water = True

            receptor_tmp = receptor.create_tmp_file(write_stream=True)
            receptor_tmp_filled = receptor.create_tmp_file()

        with capture_output() as (out, err):
            mdl = complete_pdb(env, receptor_tmp, transfer_res_num=transfer_res_num)
            mdl.write(
                file=receptor_tmp_filled,
                model_format=FileFormatHandler.get_file_ext_full_name(
                    receptor.output_fmt
                ),
            )
        if out.getvalue():
            logging.warning(out.getvalue())
        if err.getvalue():
            logging.error(err.getvalue())

        receptor.set_file(receptor_tmp_filled)
        os.remove(receptor_tmp)
        os.remove(receptor_tmp_filled)


class AddMissingAtomsOperation:
    """Add missing atoms (heavy atoms and hydrogens) to the PDB file."""

    def run_modeller(self, receptor: Receptor, transfer_res_num: bool = False) -> None:
        complete_pdb = CompletePDBOperation()
        complete_pdb.run_modeller(receptor, transfer_res_num=transfer_res_num)


class ReplaceNonStdResiduesOperation:
    """Replace non-standard residues with standard ones."""

    def change_hetatm_to_atom(
        self, receptor: Receptor, modified_res: dict[str, str]
    ) -> Receptor:
        """Change HETATM to ATOM for non-standard residues; required for MODELLER."""
        f_out = io.StringIO()
        receptor.current_file_stream.seek(0)
        for line in receptor.current_file_stream:
            resname = line[17:21].strip()
            if not line.startswith(("ATOM", "HETATM")) or resname not in modified_res:
                f_out.write(line)
            else:
                new_line = "ATOM  " + line[6:]
                f_out.write(new_line)

        receptor.current_file_stream.close()
        f_out.seek(0)
        receptor.current_file_stream = f_out
        return receptor

    def change_and_prune_non_std_residues(
        self, receptor: Receptor, nstds_to_std: dict[str, str]
    ) -> Receptor:
        """Prune non-backbone atoms from non-standard residues."""
        parser = receptor.get_biopython_parser()
        file_id = os.path.splitext(os.path.basename(receptor.file))[0]
        receptor.current_file_stream.seek(0)

        with warnings.catch_warnings(record=True) as warns:
            warnings.simplefilter("always")
            structure = parser.get_structure(file_id, receptor.current_file_stream)

        for residue in structure.get_residues():
            if residue.resname in nstds_to_std:
                residue.resname = nstds_to_std[residue.resname]  # change to standard
                for atom in residue.child_list.copy():
                    if atom.name not in ["N", "CA", "C", "O"]:
                        residue.detach_child(atom.id)

        for warn in warns:
            if warn.category == PDBExceptions.PDBConstructionWarning:
                logging.warning(f"{self.file}: {str(warn.message).split("\n")[0]}")

        receptor.close_file_stream()  # close the original file stream
        receptor.current_file_stream = io.StringIO()

        file_io = receptor.get_biopython_file_io(receptor.file_ext.strip("."))
        file_io.set_structure(structure)
        file_io.save(receptor.current_file_stream)
        return receptor

    def replace_non_std_residues(
        self, receptor: Receptor, transfer_res_num: bool = False
    ) -> None:
        complete_pdb = CompletePDBOperation()
        receptor = self.change_hetatm_to_atom(receptor, nonstd_residues.nstds_to_std)
        receptor = self.change_and_prune_non_std_residues(
            receptor, nonstd_residues.nstds_to_std
        )
        complete_pdb.run_modeller(
            receptor, transfer_res_num=transfer_res_num
        )  # this also reconstructs the missing atoms

    def run_modeller(self, receptor: Receptor, transfer_res_num: bool = False) -> None:
        self.replace_non_std_residues(receptor, transfer_res_num)

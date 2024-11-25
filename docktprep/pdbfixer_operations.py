import logging
from typing import Protocol

from openmm.app import PDBFile
from pdbfixer import PDBFixer

__all__ = ["PDBFixerOperation", "ReplaceNonStdResidues"]


class PDBFixerOperation(Protocol):
    def fix(self, fixer: PDBFixer) -> PDBFixer: ...


class ReplaceNonStdResidues:
    def fix(self, fixer: PDBFixer):
        fixer.findNonstandardResidues()
        if not fixer.nonstandardResidues:
            logging.info("No non-standard residues found.")
            return fixer

        for nst_resi, replace_resi in fixer.nonstandardResidues:
            logging.info(
                f"Replacing non-standard residue {nst_resi.name} with {replace_resi}"
            )
        fixer.replaceNonstandardResidues()
        return fixer

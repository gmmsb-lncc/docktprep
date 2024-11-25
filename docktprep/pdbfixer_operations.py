import logging
from typing import Protocol

from pdbfixer import PDBFixer

__all__ = [
    "PDBFixerOperation",
    "ReplaceNonStdResidues",
    "AddMissingHeavyAtoms",
    "AddMissingResidues",
    "AddMissingHydrogens",
]


class PDBFixerOperation(Protocol):
    def fix(self, fixer: PDBFixer) -> PDBFixer: ...


class ReplaceNonStdResidues:
    """Replace non-standard residues with their standard counterparts."""

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

        # we may need to add missing atoms for the new residues (according to the docs),
        # but the original implmentation does not check for inexistent dicts that
        # are apparently created only after calling certain methods
        fixer.missingResidues = {}
        fixer.missingAtoms = {}
        fixer.missingTerminals = {}
        fixer.addMissingAtoms()
        return fixer


class AddMissingHeavyAtoms:
    """Adds missing heavy atoms to the structure.

    There may be missing heavy atoms, e.g. in flexible regions that could not be clearly
    resolved from the electron density. Many PDB files are also missing terminal atoms
    that should be present at the ends of chains.
    """

    def fix(self, fixer: PDBFixer):
        # create attributes that are necessary for the method to work, but are created
        # elsewhere for some reason
        fixer.missingResidues = {}

        fixer.findMissingAtoms()
        if not fixer.missingAtoms and not fixer.missingTerminals:
            logging.info("No missing heavy atoms found.")
            return fixer

        for residue, missing_atoms in fixer.missingAtoms.items():
            atoms = [f"{atom.name}" for atom in missing_atoms]
            resi = f"{residue.name}{residue.id} in chain {residue.chain.id}"
            logging.info(
                f"Adding missing heavy atoms {", ".join(atoms)} to residue {resi}"
            )
        for residue, missing_atoms in fixer.missingTerminals.items():
            atoms = [f"{atom.name}" for atom in missing_atoms]
            resi = f"{residue.name}{residue.id} in chain {residue.chain.id}"
            logging.info(
                f"Adding missing heavy atoms {", ".join(atoms)} to terminal residue {resi}"
            )

        fixer.addMissingAtoms()
        return fixer


class AddMissingResidues:
    """Adds missing residues to the structure based on the SEQRES records."""

    def fix(self, fixer: PDBFixer):
        fixer.findMissingResidues()
        if not fixer.missingResidues:
            logging.info("No missing residues or no SEQRES records found.")
            return fixer

        for _, resi in fixer.missingResidues.items():
            missing_resi = ", ".join(resi)
            logging.info(f"Adding missing residues: {missing_resi}")

        # add necessary missing attributes
        fixer.missingAtoms = {}
        fixer.missingTerminals = {}
        fixer.addMissingAtoms()
        return fixer


class AddMissingHydrogens:
    """Adds missing hydrogen atoms to the structure."""

    def __init__(self, ph: float = 7.0):
        self.ph = ph

    def fix(self, fixer: PDBFixer):
        logging.info(f"Adding missing hydrogen atoms using pH {self.ph}.")
        fixer.addMissingHydrogens(self.ph)
        return fixer

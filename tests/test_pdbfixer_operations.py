import pytest
from pdbfixer import PDBFixer

from docktprep.pdbfixer_operations import *


def test_no_std_residues_replaces():
    fixer = PDBFixer("tests/data/1bkx.pdb")
    fixer.findNonstandardResidues()
    assert fixer.nonstandardResidues

    replace = ReplaceNonStdResidues()
    fixer = replace.fix(fixer)
    fixer.findNonstandardResidues()
    assert not fixer.nonstandardResidues


def test_no_std_residues_does_not_raise():
    fixer = PDBFixer("tests/data/1az5.pdb")
    replace = ReplaceNonStdResidues()
    fixer = replace.fix(fixer)
    assert fixer


def test_adds_missing_heavy_atoms():
    fixer = PDBFixer("tests/data/1az5.pdb")
    fixer.missingResidues = []
    fixer.findMissingAtoms()
    assert fixer.missingAtoms
    # assert fixer.missingTerminals

    add_missing = AddMissingHeavyAtoms()
    fixer = add_missing.fix(fixer)
    fixer.findMissingAtoms()
    assert not fixer.missingAtoms
    assert not fixer.missingTerminals


def test_missing_heavy_atoms_does_not_raise():
    fixer = PDBFixer("tests/data/9ins.pdb")
    add_missing_heavy = AddMissingHeavyAtoms()
    fixer = add_missing_heavy.fix(fixer)
    assert fixer


def test_adds_missing_residues():
    fixer = PDBFixer("tests/data/1az5.pdb")
    fixer.missingResidues = {}
    fixer.findMissingResidues()
    assert fixer.missingResidues

    add_missing = AddMissingResidues()
    fixer = add_missing.fix(fixer)
    fixer.findMissingResidues()
    assert not fixer.missingResidues


def test_missing_residues_does_not_raise():
    fixer = PDBFixer("tests/data/9ins.pdb")
    add_missing_residues = AddMissingResidues()
    fixer = add_missing_residues.fix(fixer)
    assert fixer


def test_add_missing_hydrogens_does_not_raise():
    fixer = PDBFixer("tests/data/1az5.pdb")
    add_missing_hs = AddMissingHydrogens()
    fixer = add_missing_hs.fix(fixer)
    assert fixer

import pytest
from pdbfixer import PDBFixer

from docktprep.pdbfixer_operations import *


def test_no_std_residues_replace():
    fixer = PDBFixer("tests/data/1bkx.pdb")
    fixer.findNonstandardResidues()
    assert fixer.nonstandardResidues

    replace = ReplaceNonStdResidues()
    fixer = replace.fix(fixer)
    fixer.findNonstandardResidues()
    assert not fixer.nonstandardResidues

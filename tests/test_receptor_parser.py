import io
import tempfile
from copy import copy

import pytest

from docktprep.pdbfixer_operations import *
from docktprep.receptor_parser import PDBSanitizerFactory, Receptor


def test_pdb_receptor_open_file_stream():
    receptor = Receptor("tests/data/1az5.pdb")
    assert receptor.current_file_stream is not None
    assert type(receptor.current_file_stream) == io.TextIOWrapper
    receptor.close_file_stream()


def test_pdb_receptor_open_file_stream_error():
    filename = "notafile.pdb"
    with pytest.raises(FileNotFoundError):
        Receptor(filename)


def test_pdb_receptor_close_file_stream():
    receptor = Receptor("tests/data/1az5.pdb")
    receptor.close_file_stream()
    with pytest.raises(ValueError):
        receptor.current_file_stream.read()


def test_sanitize_file():
    receptor = Receptor("tests/data/1az5.pdb")
    receptor.sanitize_file()
    assert receptor.current_file_stream is not None
    assert type(receptor.current_file_stream) == io.StringIO
    receptor.close_file_stream()


def test_sanitize_file_model_id_does_not_exist():
    sanitizer_factory = PDBSanitizerFactory(model_id=999)
    receptor = Receptor("tests/data/1az5.pdb", sanitizer=sanitizer_factory)

    with pytest.raises(ValueError):
        receptor.sanitize_file()
    receptor.close_file_stream()


def test_write_file_stream():
    receptor = Receptor("tests/data/9ins.pdb")
    receptor.sanitize_file()
    with tempfile.NamedTemporaryFile(delete=True) as tmp:
        current_file_stream = copy(receptor.current_file_stream)
        receptor.write_and_close_file_stream(tmp.name)
        with open(tmp.name, "r") as f:
            current_file_stream.seek(0)
            assert f.read() == current_file_stream.read()

    receptor.close_file_stream()


def test_sanitize_and_fix_structure_does_not_raise():
    receptor = Receptor("tests/data/1az5.pdb")
    fix_ops = [
        ReplaceNonStdResidues(),
        AddMissingHeavyAtoms(),
        AddMissingResidues(),
        AddMissingHydrogens(),
    ]
    receptor.sanitize_file()
    receptor.fix_structure(fix_ops)
    with tempfile.NamedTemporaryFile(delete=True) as tmp:
        receptor.write_and_close_file_stream(tmp.name)
        with open(tmp.name, "r") as f:
            assert f.read()
    receptor.close_file_stream()


# def test_sanitize_file_adds_seqres_to_stream():
#     receptor = Receptor("tests/data/1az5.pdb")
#     seqres = receptor.get_seqres_from_stream()
#     assert seqres
#     receptor.current_file_stream.seek(0)
#     receptor.sanitize_file()  # should add SEQRES records to the top of the stream
#     receptor.current_file_stream.seek(0)
#     assert "SEQRES" in receptor.current_file_stream.getvalue()

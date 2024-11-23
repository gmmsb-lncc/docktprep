import io
import tempfile
from unittest.mock import Mock

import pytest

from docktprep.receptor_parser import PDBSanitizer, PDBSanitizerFactory, Receptor


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
        receptor.write_file_stream(tmp.name)
        with open(tmp.name, "r") as f:
            receptor.current_file_stream.seek(0)
            assert f.read() == receptor.current_file_stream.read()

    receptor.close_file_stream()

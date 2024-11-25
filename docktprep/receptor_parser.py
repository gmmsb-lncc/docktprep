import io
import logging
import os
import warnings

from Bio.PDB import PDBIO, PDBExceptions, PDBParser, Structure
from Bio.PDB.PDBIO import Select
from openmm.app import PDBFile
from pdbfixer import PDBFixer

from .pdbfixer_operations import PDBFixerOperation


class PDBSanitizer(Select):
    def __init__(self, model_id) -> None:
        """Select objects from a PDB structure before writing.

        Parameters
        ----------
        model_id : int
            Model ID to select, by default 0
        """
        super().__init__()
        self.model_id = model_id

    def setup_structure(self, structure: Structure):
        """Setup the structure and collect disordered atoms."""
        self.structure = structure
        self.structure_model_ids = [m.id for m in structure.get_models()]
        self.disordered_atoms = []

        # collect disordered atoms
        for atom in structure.get_atoms():
            if atom.is_disordered():
                self.disordered_atoms.append(atom)

    def disordered_atoms_to_reject_by_occupancy(self):
        """Return a list of disordered atom ids (number) to reject based on occupancy."""
        reject_ids = []
        for atom in self.disordered_atoms:
            highest_occ_altloc = atom.get_altloc()  # defaults to first altloc
            for child in atom.child_dict.values():
                if child.get_altloc() != highest_occ_altloc:
                    reject_ids.append(child.get_serial_number())
        return reject_ids

    def accept_model(self, model):
        if not self.model_id in self.structure_model_ids:
            e = f"Model ID {self.model_id} not found in structure."
            logging.error(e)
            raise ValueError(e)

        if model.id == self.model_id:
            logging.info(f"Selected model ID {self.model_id}")
            return True
        return False

    def accept_atom(self, atom):
        """If atom is disordered, select the highest occupancy atom."""
        reject_ids = self.disordered_atoms_to_reject_by_occupancy()
        if atom.get_serial_number() in reject_ids:
            logging.info(
                f"Ignoring lower occupancy atom ({atom.get_serial_number()} {atom.get_name()})"
            )
            return False
        return True


class PDBSanitizerFactory:
    def __init__(self, model_id: int = 0):
        self.model_id = model_id

    def create_sanitizer(self, structure: Structure) -> PDBSanitizer:
        sanitizer = PDBSanitizer(model_id=self.model_id)
        sanitizer.setup_structure(structure)
        return sanitizer


class Receptor:
    def __init__(
        self,
        file: str,
        output_fmt: str = "",
        sanitizer: PDBSanitizerFactory = PDBSanitizerFactory(),  # sanitizer with default args
    ) -> None:
        self.file = file
        self.file_ext = os.path.splitext(file)[1]
        self.sanitizer = sanitizer
        self.current_file_stream = self.open_file_stream()
        self.output_fmt = output_fmt if output_fmt else self.file_ext.strip(".")

    def open_file_stream(self) -> io.TextIOWrapper:
        try:
            return open(self.file, "r")
        except FileNotFoundError as e:
            logging.error(f"File not found: {self.file}")
            raise e

    def close_file_stream(self):
        self.current_file_stream.close()

    def write_and_close_file_stream(self, file: str):
        self.current_file_stream.seek(0)
        with open(file, "w") as f:
            f.write(self.current_file_stream.read())
        self.close_file_stream()

    def get_seqres_from_stream(self):
        self.current_file_stream.seek(0)
        seqres = []
        for line in self.current_file_stream:
            if line.startswith("SEQRES"):
                seqres.append(line)
        return seqres

    def add_seqres_to_stream(self, seqres: list[str]):
        """Adds SEQRES records to the beginning of the current file stream."""
        seqres_str = "".join(seqres)

        self.current_file_stream.seek(0)
        original_content = self.current_file_stream.read()

        self.current_file_stream = io.StringIO()
        self.current_file_stream.write(seqres_str + original_content)

    def get_biopython_parser(self):
        if self.file_ext == ".pdb":
            return PDBParser(PERMISSIVE=True, QUIET=False)
        else:
            e = f"Unsupported file extension: {self.file_ext}"
            logging.error(e)
            raise ValueError(e)

    def get_pdbfixer_parser(self):
        if self.file_ext == ".pdb":
            return PDBFixer
        else:
            e = f"Unsupported file extension: {self.file_ext}"
            logging.error(e)
            raise ValueError(e)

    def get_biopython_file_io(self, output_fmt: str = ""):
        if output_fmt == "pdb":
            return PDBIO()
        else:
            e = f"Unsupported output format: {self.file_ext}"
            logging.error(e)
            raise ValueError(e)

    def sanitize_file(self) -> None:
        """Sanitize the receptor file using biopython.

        Catches common PDB exceptions and errors. Save the sanitized file
        to a new `current_file_stream`. Logs any warnings.
        """
        seqres = self.get_seqres_from_stream()
        parser = self.get_biopython_parser()
        file_id = os.path.splitext(os.path.basename(self.file))[0]
        self.current_file_stream.seek(0)

        with warnings.catch_warnings(record=True) as warns:
            warnings.simplefilter("always")
            structure = parser.get_structure(file_id, self.current_file_stream)

        for warn in warns:
            if warn.category == PDBExceptions.PDBConstructionWarning:
                logging.warning(f"{self.file}: {str(warn.message).split("\n")[0]}")

        self.close_file_stream()  # close the original file stream
        self.current_file_stream = io.StringIO()

        file_io = self.get_biopython_file_io(self.file_ext.strip("."))
        file_io.set_structure(structure)
        sanitizer = self.sanitizer.create_sanitizer(structure)
        file_io.save(self.current_file_stream, write_end=True, select=sanitizer)
        # add SEQRES records to the beginning of the file, since bioptyhon excludes
        # the header info
        # self.add_seqres_to_stream(seqres)

    def fix_structure(self, fix_ops: list[PDBFixerOperation]) -> None:
        self.current_file_stream.seek(0)
        self.current_file_stream.seek(0)
        fixer_parser = self.get_pdbfixer_parser()
        fixer = fixer_parser(pdbfile=self.current_file_stream)
        for op in fix_ops:
            fixer = op.fix(fixer)

        self.close_file_stream()  # close the original file stream
        self.current_file_stream = io.StringIO()
        PDBFile.writeFile(fixer.topology, fixer.positions, self.current_file_stream)

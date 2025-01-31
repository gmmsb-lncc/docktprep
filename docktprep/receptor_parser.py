import io
import logging
import os
import tempfile
import warnings

from Bio.PDB import PDBIO, PDBExceptions, PDBParser, Structure
from Bio.PDB.PDBIO import Select


class PDBSanitizer(Select):
    def __init__(
        self,
        model_id: int = 0,
        remove_disorder: bool = True,
        remove_hetresi: bool = True,
        remove_water: bool = True,
        **kwargs,
    ) -> None:
        """Biopython Select class to sanitize PDB files.

        Parameters
        ----------
        model_id : int
            model id to select (default: 0)
        remove_disorder : bool
            remove disordered atoms (default: True)
        remove_hetatms : bool
            remove HETATM atoms (default: True)
        """
        super().__init__()
        self.model_id = model_id
        self.remove_disorder = remove_disorder
        self.remove_hetresi = remove_hetresi
        self.remove_water = remove_water

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
        if atom.get_serial_number() in reject_ids and self.remove_disorder:
            logging.info(
                f"Ignoring lower occupancy atom ({atom.get_serial_number()} {atom.get_name()})"
            )
            return False

        if atom.get_parent().id[0] == "W" and self.remove_water:
            logging.info(
                f"Removing WATER atom ({atom.get_serial_number()} {atom.get_name()})"
            )
            return False

        if "H_" in atom.get_parent().id[0] and self.remove_hetresi:
            logging.info(
                f"Removing HETATM atom ({atom.get_serial_number()} {atom.get_name()})"
            )
            return False
        return True


class PDBSanitizerFactory:
    def __init__(self, **kwargs) -> None:
        self.kwargs = kwargs

    def create_sanitizer(self, structure: Structure) -> PDBSanitizer:
        sanitizer = PDBSanitizer(**self.kwargs)
        sanitizer.setup_structure(structure)  # required logic to setup the sanitizer
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

    def set_file(self, file: str):
        self.current_file_stream.close()
        self.file = file
        self.file_ext = os.path.splitext(file)[1]
        self.current_file_stream = self.open_file_stream()

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

    def get_biopython_parser(self):
        if self.file_ext == ".pdb":
            return PDBParser(PERMISSIVE=True, QUIET=False)
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

    def create_tmp_file(self, write_stream: bool = False) -> str:
        with tempfile.NamedTemporaryFile(delete=False) as tmp:
            tmp_file = tmp.name
            if write_stream:
                self.write_and_close_file_stream(tmp_file)
            return tmp_file

    def sanitize_file(self) -> None:
        """Sanitize the receptor file using biopython.

        Catches common PDB exceptions and errors. Save the sanitized file
        to a new `current_file_stream`. Logs any warnings.
        """
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

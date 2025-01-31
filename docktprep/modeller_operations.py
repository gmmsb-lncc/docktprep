import logging
import os
from typing import Protocol

from modeller import *
from modeller.scripts import complete_pdb

from .receptor_parser import Receptor
from .stdout_manager import capture_output, suppress_output


class ModellerOperation(Protocol):
    def run_modeller(self, receptor: Receptor) -> None: ...


class AddMissingAtomsOperation:
    def add_missing_atoms(self, receptor: Receptor) -> None:
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
            mdl = complete_pdb(env, receptor_tmp)
            mdl.write(file=receptor_tmp_filled)
        if out.getvalue():
            logging.warning(out.getvalue())
        if err.getvalue():
            logging.error(err.getvalue())

        receptor.set_file(receptor_tmp_filled)
        os.remove(receptor_tmp)
        os.remove(receptor_tmp_filled)

    def run_modeller(self, receptor: Receptor) -> None:
        self.add_missing_atoms(receptor)

import argparse
import os

from docktprep.pdbfixer_operations import *
from docktprep.receptor_parser import PDBSanitizerFactory, Receptor

from .logs import configure_logging


def main(args):
    sanitizer = PDBSanitizerFactory(model_id=args.select_model)
    receptor = Receptor(
        args.receptor,
        sanitizer=sanitizer,
    )

    fixer_ops = []
    if args.add_missing_residues:
        fixer_ops.append(AddMissingResidues())
    if args.replace_non_standard:
        fixer_ops.append(ReplaceNonStdResidues())
    if args.add_missing_atoms:
        fixer_ops.append(AddMissingHeavyAtoms())
    if args.add_hydrogens:
        fixer_ops.append(AddMissingHydrogens(ph=args.pH))

    receptor.sanitize_file()
    receptor.fix_structure(fixer_ops)
    receptor.write_and_close_file_stream(args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="DockTPrep: Create DockThor input files from PDB, mmCIF or MOL2 formats.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-r",
        "--receptor",
        help="Receptor file in PDB or mmCIF format.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--output",
        help="Output file name to save the prepared structure.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--log-output",
        type=str,
        help="Output file for logging.",
        default=None,
    )
    receptor_operations = parser.add_argument_group("Receptor options")
    receptor_operations.add_argument(
        "--select-model",
        metavar="MODEL_IDX",
        type=int,
        default=0,
        help="Select a model from the input file, in case of multiple models.",
    )
    # receptor_operations.add_argument(
    #     "--convert-to-pdbx",
    #     action="store_true",
    #     help="Convert receptor input file to PDBx/mmCIF format.",
    # )
    receptor_operations.add_argument(
        "--replace-non-standard",
        action="store_true",
        help="Replace non-standard residues with their standard counterparts.",
    )
    receptor_operations.add_argument(
        "--add-missing-atoms",
        action="store_true",
        help="Add missing heavy atoms to the structure.",
    )
    receptor_operations.add_argument(
        "--add-missing-residues",
        action="store_true",
        help="Add missing residues to the structure based on SEQRES records.",
    )
    receptor_operations.add_argument(
        "--add-hydrogens",
        action="store_true",
        help="Add missing hydrogen atoms to the structure using the specified pH.",
    )
    receptor_operations.add_argument(
        "--pH",
        type=float,
        default=7.0,
        help="pH value for adding missing hydrogen atoms to the receptor structure.",
    )

    args = parser.parse_args()

    #
    log_output_file = (
        f"docktprep_{os.path.basename(args.receptor)}.log"
        if args.log_output is None
        else args.log_output
    )
    configure_logging(output_file=log_output_file)
    main(args)

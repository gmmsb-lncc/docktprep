import argparse
import os

from docktprep.receptor_parser import PDBSanitizerFactory, Receptor

from .logs import configure_logging


def main():
    args = configure_argparser()
    configure_logging(args.log_file)

    # receptor parsing
    sanitizer = PDBSanitizerFactory(
        model_id=args.sel_model,
        remove_hetresi=args.remove_hetresi,
        remove_water=args.remove_water,
    )

    receptor = Receptor(
        args.receptor,
        sanitizer=sanitizer,
    )

    receptor.sanitize_file()
    receptor.write_and_close_file_stream(args.output)


def configure_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="DockTPrep: Prepare protein-ligand structures and create DockThor input files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-r",
        "--receptor",
        help="Receptor file in PDB format.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output file name to save the prepared structure.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--log-file",
        type=str,
        help="Output file for logging.",
        default=None,
    )
    receptor_operations = parser.add_argument_group("receptor options")

    receptor_operations.add_argument(
        "--remove-hetresi",
        action="store_true",
        help="Remove HETATM records.",
    )

    receptor_operations.add_argument(
        "--remove-water",
        action="store_true",
        help="Remove water molecules.",
    )

    receptor_operations.add_argument(
        "--sel-model",
        type=int,
        default=0,
        help="Select a model from the input file using its index, in case of multiple models.",
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()

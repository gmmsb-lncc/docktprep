import argparse
import os

from docktprep.receptor_parser import PDBSanitizerFactory, Receptor

from .logs import configure_logging


def main():
    args = configure_argparser()
    configure_logging(args.log_output)
    sanitizer = PDBSanitizerFactory(model_id=args.select_model)
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
    args = parser.parse_args()

    #
    log_output_file = (
        f"docktprep_{os.path.basename(args.receptor)}.log"
        if args.log_output is None
        else args.log_output
    )
    args.log_output = log_output_file
    return args


if __name__ == "__main__":
    main()

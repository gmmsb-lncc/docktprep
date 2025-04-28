import argparse

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

    # modeller operations
    receptor = modeller_operations(receptor, args)

    # write receptor to output file
    receptor.write_and_close_file_stream(args.output)


def modeller_operations(receptor: Receptor, args: argparse.Namespace):
    try:
        from docktprep import modeller_operations
    except ImportError:
        raise ImportError(f"MODELLER is required to use this feature.")

    mdlops = list()
    if args.add_missing_atoms and args.replace_nstd_res:
        # this will also add missing atoms
        mdlops.append(modeller_operations.ReplaceNonStdResiduesOperation())
    elif args.add_missing_atoms:
        mdlops.append(modeller_operations.AddMissingAtomsOperation())
    elif args.replace_nstd_res:
        mdlops.append(modeller_operations.ReplaceNonStdResiduesOperation())

    for mdlop in mdlops:
        mdlop.run_modeller(receptor, transfer_res_num=args.transfer_res_num)

    return receptor


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

    # receptor operations
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
        help="Select a model from the input file using its index.",
    )
    receptor_operations.add_argument(
        "--add-missing-atoms",
        action="store_true",
        help="Add missing heavy and hydrogen atoms (requires MODELLER).",
    )
    receptor_operations.add_argument(
        "--replace-nstd-res",
        action="store_true",
        help="Replace non-standard residues with their standard counterparts (requires MODELLER).",
    )
    receptor_operations.add_argument(
        "--transfer-res-num",
        action="store_true",
        help="Retains the residue numbering from the original PDB (MODELLER).",
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()

"""Configure application logging."""

import logging


def configure_logging(
    output_file: str = "docktprep.log", level: int = logging.INFO
) -> None:
    """Configure application logging."""
    logging.basicConfig(
        filename=output_file,
        filemode="a",
        level=level,
        format="%(asctime)s: %(levelname)s: %(message)s",
    )

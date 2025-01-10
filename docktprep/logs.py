"""Configure application logging."""

import logging


def configure_logging(
    output_file: str | None = None, level: int = logging.INFO
) -> None:
    """Configure application logging."""
    logging.basicConfig(
        filename=output_file if output_file else None,
        filemode="w",
        level=level,
        format="%(asctime)s: %(levelname)s: %(message)s",
    )

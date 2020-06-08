"""
"""
import logging
from pathlib import Path
from importlib_resources import open_text

logger = logging.getLogger(__name__)

CONFIGURATION_FILE_NAME = "naibr.config"


def add_arguments(parser):
    parser.add_argument("analysis_folder", type=Path, help="BLR analysis folder.")
    parser.add_argument("naibr_out", type=Path, help="Config: NAIBR output folder.")
    parser.add_argument("bamfile", type=Path, help="Config: Bam file for which LSVs should be called.")
    parser.add_argument("--threads", type=int, default=1, help="Config: Number of threads for which to run NAIBR.")


def main(args):

    new_value = {"bam_file": str(args.bamfile.resolve()), "outdir": str(args.naibr_out.resolve()), "threads": str(args.threads)}
    mod_keys = tuple(new_value.keys())
    configuration = open_text("blr", CONFIGURATION_FILE_NAME)
    with (args.analysis_folder / CONFIGURATION_FILE_NAME).open("w") as f:
        for row in configuration:
            if row.startswith(mod_keys):
                key = row.split("=", maxsplit=1)[0]
                row = key + "=" + new_value[key] + "\n"
                logger.info(f"Setting config value '{row.strip()}'")
            f.write(row)

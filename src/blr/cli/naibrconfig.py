"""
Builds NAIBR config files.

File format:
 - comments: # my documentation
 - blank rows allowed
 - variable=value
"""
import logging
from pathlib import Path
from importlib_resources import open_text

logger = logging.getLogger(__name__)

CONFIGURATION_FILE_NAME = "naibr.config"


def add_arguments(parser):
    parser.add_argument("analysis_folder", type=Path, help="BLR analysis folder.")
    parser.add_argument("bamfile", type=Path, help="Config: Bam file for which LSVs should be called.")
    parser.add_argument("naibr_out", type=Path, help="Config: NAIBR output folder.")
    parser.add_argument("--threads", type=int, default=1, help="Config: Number of threads for which to run NAIBR.")


def main(args):
    new_value = {"bam_file": str(args.bamfile), "outdir": str(args.naibr_out),
                 "threads": str(args.threads)}
    mod_keys = tuple(new_value.keys())
    copy_and_mod_config(args.analysis_folder, CONFIGURATION_FILE_NAME, new_value, mod_keys)


def copy_and_mod_config(analysis_folder: Path, file_name: Path, new_value: dict, mod_keys: tuple):
    configuration = open_text("blr", file_name)
    with (analysis_folder / file_name).open("w") as f:
        for row in configuration:
            if row.startswith(mod_keys):
                row = change_row(row, new_value)
            f.write(row)


def change_row(row, new_value):
    key = row.split("=", maxsplit=1)[0]
    row = key + "=" + new_value[key] + "\n"
    logger.info(f"Setting config value '{row.strip()}'")
    return row

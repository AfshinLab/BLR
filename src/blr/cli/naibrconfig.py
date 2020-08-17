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
import os

logger = logging.getLogger(__name__)

CONFIGURATION_FILE_NAME = "naibr.config"


def add_arguments(parser):
    parser.add_argument(
        "--bam-file", required=True, type=Path,
        help="Input BAM for ananlysis.")
    parser.add_argument(
        "--outdir", default=Path(os.getcwd()), type=Path,
        help="NAIBR output directory. Default: %(default)s")
    parser.add_argument(
        "-d", "--distance", default=10000, type=int,
        help="Maximum distance between read-pairs in a linked-read. Default: %(default)s")
    parser.add_argument(
        "--min-mapq", default=40, type=int,
        help="Minimum mapping quality. Default: %(default)s")
    parser.add_argument(
        "--min-sv", type=int,
        help="Minimum size of structural variant. Default: estimated from input bam")
    parser.add_argument(
        "--threads", type=int, default=1,
        help="Number of threads for which to run NAIBR. Default: %(default)s")
    parser.add_argument(
        "-k", "--min-overlaps", type=int, default=3,
        help="Minimum number of barcode overlaps supporting a candidate NA. Default: %(default)s")
    parser.add_argument(
        "-o", "--output", type=Path, required=True,
        help="Name of output config file.")


def main(args):
    parameters = {
        "min_mapq": args.min_mapq,
        "bam_file": args.bam_file.resolve(),
        "outdir": args.outdir.resolve(),
        "d": args.distance,
        "min_sv": args.min_sv,
        "threads": args.threads,
        "k": args.min_overlaps
    }

    copy_and_mod_config(CONFIGURATION_FILE_NAME, args.output, parameters)


def copy_and_mod_config(template_file: str, output_file: Path, parameters):
    configuration = open_text("blr", template_file)
    with output_file.open("w") as f:
        for row in configuration:
            if not row.startswith("#") and "=" in row:
                row = change_row(row, parameters)
            f.write(row)


def change_row(row, parameters):
    key = row.split("=", maxsplit=1)[0]
    row = key + "=" + str(parameters[key]) + "\n"
    logger.info(f"Setting config value '{row.strip()}'")
    return row

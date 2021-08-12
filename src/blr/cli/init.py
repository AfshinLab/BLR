"""
Create and initialize a new analysis directory.
"""
from importlib_resources import read_binary
import logging
import os
import os.path
from pathlib import Path
import sys
from typing import List

from blr.utils import guess_paired_path, ACCEPTED_LIBRARY_TYPES
from blr.cli.config import change_config, load_yaml

logger = logging.getLogger(__name__)


CONFIGURATION_FILE_NAME = "blr.yaml"
MULTIQC_CONFIG_FILE_NAME = "multiqc_config.yaml"
KEY_FILES = {"final.bam", "final.molecule_stats.filtered.tsv", "barcodes.clstr.gz"}  # For dbs (blr) and tellseq libs
KEY_FILES2 = {"final.bam", "final.molecule_stats.filtered.tsv"}  # For 10x and stLFR libs


def add_arguments(parser):
    required = parser.add_mutually_exclusive_group()
    required.add_argument(
        "--reads1",
        "--r1",
        type=Path,
        metavar="READS",
        help="First paired-end read file (.fastq.gz). The second is found automatically.",
    )
    parser.add_argument(
        "-l",
        "--library-type",
        required=True,
        choices=ACCEPTED_LIBRARY_TYPES,
        help="Select library type from currently available technologies: %(choices)s."
    )
    required.add_argument(
        "-w",
        "--from-workdir",
        type=Path,
        metavar="DIR",
        action="append",
        help="Initailize new analysis directory based on previous analysis instead of FASTQ reads. Will identify key "
             "files from previous run(s) and use the for the basis of new analysis. If multiple DIRs are provided "
             "the files will be merged appropriately."
    )
    parser.add_argument("directory", type=Path, help="New analysis directory to create")


def main(args):
    # TODO remove in next version
    if args.library_type == "blr":
        logger.warn("Deprecation: Replacing 'blr' with 'dbs' for library_type.")
        args.library_type = "dbs"

    if args.reads1:
        init(args.directory, args.reads1, args.library_type)
    elif args.from_workdir:
        init_from_dir(args.directory, args.from_workdir, args.library_type)


def init(directory: Path, reads1: Path, library_type: str):
    if " " in str(directory):
        logger.error("The name of the analysis directory must not contain spaces")
        sys.exit(1)

    fail_if_inaccessible(reads1)
    reads2 = guess_paired_path(reads1)
    if reads2 is None:
        logger.error("Could not determine second file of paired-end reads")
        sys.exit(1)
    fail_if_inaccessible(reads2)

    create_and_populate_analysis_directory(directory, reads1, reads2, library_type)

    logger.info(f"Directory {directory} initialized.")
    logger.info(
        'Edit %s/%s, then run "cd %s && blr run" to start the analysis',
        directory,
        CONFIGURATION_FILE_NAME,
        directory,
    )


def create_and_populate_analysis_directory(directory: Path, reads1: Path, reads2: Path, library_type: str):
    try_mkdir(directory)

    # Write the configuration files
    write_config_to_dir(CONFIGURATION_FILE_NAME, directory)
    write_config_to_dir(MULTIQC_CONFIG_FILE_NAME, directory)

    # Update with library type into
    change_config(directory / CONFIGURATION_FILE_NAME, [("library_type", library_type)])

    create_symlink(reads1, directory, "reads.1.fastq.gz")
    create_symlink(reads2, directory, "reads.2.fastq.gz")


def write_config_to_dir(file_name: str, directory: Path):
    # Write the configuration file
    configuration = read_binary("blr", file_name)
    with (directory / file_name).open("wb") as f:
        f.write(configuration)


def fail_if_inaccessible(path):
    try:
        with path.open():
            pass
    except OSError as e:
        logger.error("Could not open %r: %s", path, e)
        sys.exit(1)


def create_symlink(readspath, dirname, target):
    if not os.path.isabs(readspath):
        src = os.path.relpath(readspath, dirname)
    else:
        src = readspath
    os.symlink(src, os.path.join(dirname, target))


def try_mkdir(directory: Path):
    try:
        directory.mkdir()
    except OSError as e:
        logger.error(e)
        sys.exit(1)


def init_from_dir(directory: Path, workdirs: List[Path], library_type: str):
    if not all(w.is_dir() and (w / CONFIGURATION_FILE_NAME).exists() for w in workdirs):
        logger.error(f"The workdir paths must lead to directories and contain the file '{CONFIGURATION_FILE_NAME}'")
        sys.exit(1)

    required_files = {
        "blr": KEY_FILES,
        "dbs": KEY_FILES,
        "tellseq": KEY_FILES,
        "stlfr": KEY_FILES2,
        "10x": KEY_FILES2
        }[library_type]

    for file in required_files:
        if not all((w / file).exists() for w in workdirs):
            logger.error(f"The workdirs must contain the file '{file}'")
            sys.exit(1)

    configs = list(get_configs(workdirs))
    if len({c["sample_nr"] for c in configs}) != len(workdirs):
        logger.warning("The sample_nr should be different for each analysis run in order to not merge unrelated"
                       "barcodes")

    # TODO Enable re-tagging files is share same 'sample_nr'?

    library_types = [c["library_type"].replace("blr", "dbs") for c in configs]  # TODO remove replace
    if any(lib == library_type for lib in library_types):
        sys.exit(f"Trying to merge libraries of different types or not matching requested type '{library_type}':"
                 f" {', '.join([f'{dir} = {lib}' for dir, lib in zip(workdirs, library_types)])}")

    if " " in str(directory):
        logger.error("The name of the analysis directory must not contain spaces")
        sys.exit(1)

    try_mkdir(directory)

    # Write the configuration files
    write_config_to_dir(CONFIGURATION_FILE_NAME, directory)
    write_config_to_dir(MULTIQC_CONFIG_FILE_NAME, directory)

    # Update with library type into
    change_config(directory / CONFIGURATION_FILE_NAME, [("library_type", library_type)])

    input_dir = directory / "inputs"
    input_dir.mkdir()

    for nr, workdir in enumerate(workdirs, start=1):
        for file in required_files:
            target_name = f"dir{nr}.{file}"
            create_symlink(workdir / file, input_dir, target_name)

    logger.info(f"Directory {directory} initialized.")
    logger.info(f"Edit {directory}/{CONFIGURATION_FILE_NAME}.")
    logger.info(f"Run 'cd {directory} && blr run anew' to start the analysis.")


def get_configs(workdirs: List[Path]):
    for w in workdirs:
        config_path = w / CONFIGURATION_FILE_NAME
        configs, _ = load_yaml(config_path)
        yield configs

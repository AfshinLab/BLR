"""
Update configuration file. If no --set option is given the current settings are printed.
"""
import sys
import logging
from ruamel.yaml import YAML
from snakemake.utils import validate
from importlib_resources import path as resource_path
from pathlib import Path
from typing import List, Tuple
from shutil import get_terminal_size

logger = logging.getLogger(__name__)
DEFAULT_PATH = Path("blr.yaml")
SCHEMA_FILE = "config.schema.yaml"


# Script is based on repos NBISSweden/IgDisover config script.
# Link https://github.com/NBISweden/IgDiscover/blob/master/src/igdiscover/cli/config.py

def main(args):
    run(yaml_file=args.file, changes_set=args.set)


def run(yaml_file=DEFAULT_PATH, changes_set=None):
    if changes_set:
        change_config(yaml_file, changes_set)
    else:
        print_config(yaml_file)


def print_config(filename: Path):
    """
    Print out current configs to terminal.
    """
    configs, yaml = load_yaml(filename)
    width, _ = get_terminal_size()
    header = f" CONFIGS IN: {filename} "
    padding = int((width - len(header)) / 2) * "="

    # Print out current setting
    print(f"{padding}{header}{padding}")
    yaml.dump(configs, stream=sys.stdout)
    print(f"{'=' * width}")


def change_config(filename: Path, changes_set: List[Tuple[str, str]]):
    """
    Change config YAML file at filename using the changes_set key-value pairs.
    :param filename: Path to YAML config file to change.
    :param changes_set: changes to incorporate.
    """
    # Get configs from file.
    configs, yaml = load_yaml(filename)

    # Update configs
    for key, value in changes_set:
        # Convert relative paths to absolute
        value = make_paths_absolute(value, workdir=filename.parent)
        value = YAML(typ='safe').load(value)
        logger.info(f"Changing value of '{key}': {configs[key]} --> {value}.")
        item = configs

        # allow nested keys
        keys = key.split('.')
        for i in keys[:-1]:
            item = item[i]
        item[keys[-1]] = value

    # Confirm that configs is valid.
    with resource_path('blr', SCHEMA_FILE) as schema_path:
        validate(configs, str(schema_path))

    # Write first to temporary file then overwrite filename.
    tmpfile = Path(str(filename) + ".tmp")
    with open(tmpfile, "w") as file:
        yaml.dump(configs, stream=file)
    tmpfile.rename(filename)


def load_yaml(filename: Path):
    """
    Load YAML file and return the yaml object and data.
    :param filename: Path to YAML file
    :return: (data, yaml).
    """
    with open(filename) as file:
        yaml = YAML()
        data = yaml.load(file)
    return data, yaml


def make_paths_absolute(value: str, workdir: Path = Path.cwd()) -> str:
    """
    Detect if value is a relative path and make it absolut if so.
    :param value: Parameter value from arguments
    :param workdir: Path to workdir. Default: CWD
    :return:
    """
    if "../" in value and (workdir / value).exists():
        if (workdir / value).is_symlink():
            return str((workdir / value).absolute())
        return str((workdir / value).resolve(True))
    return value


def add_arguments(parser):
    parser.add_argument("-s", "--set", nargs=2, metavar=("KEY", "VALUE"), action="append",
                        help="Set KEY to VALUE. Use KEY.SUBKEY[.SUBSUBKEY...] for nested keys. For empty values "
                             "write 'null'. Can be given multiple times.")
    parser.add_argument("-f", "--file", default=DEFAULT_PATH, type=Path,
                        help="Configuration file to modify. Default: %(default)s in current directory.")

"""
Run the BLR pipeline.

This is a small wrapper around Snakemake that sets some default parameters.

To run full pipeline use:

    $ blr run

To do a dry-run (not executing anything, only showing what would be done) use:

    $ blr run -n

For additional arguments related to Snakemake run '$ snakemake -h' or look at the official
documentation at https://snakemake.readthedocs.io/en/stable/executing/cli.html.
"""

# Snakemake wrapping parially based on:
#  - http://ivory.idyll.org/blog/2020-improved-workflows-as-applications.html
#  - https://github.com/dib-lab/charcoal/blob/latest/charcoal/__main__.py

from importlib_resources import path as resource_path
import logging
import sys
import os
import subprocess
from typing import List

from snakemake.utils import available_cpu_count

logger = logging.getLogger(__name__)


def add_arguments(parser):
    parser.add_argument(
        '-c', '--cores', metavar='N', type=int, default=available_cpu_count(),
        help='Run on at most N CPU cores in parallel. '
        'Default: %(default)s (all available cores).'
    )
    parser.add_argument(
        '--anew', action="store_true", default=False,
        help="Use if initializing from previous analysis run(s)."
    )
    parser.add_argument(
        '--no-use-conda', action="store_true", default=False,
        help="Skip passing argument '--use-conda' to snakemake."
    )

    # This argument will not capture **all** additional arguments. Instead parse_known_args()
    # is used in __main__.py to add any arguments not captured here to snakemake_args.
    smk_args = parser.add_argument_group("snakemake arguments")
    smk_args.add_argument(
        'snakemake_args', nargs="*",
        help="Arguments passed to snakemake. For info about snakemake options run "
             "'snakemake --help'."
    )


def main(args):
    try:
        if args.anew:
            run(cores=args.cores,
                no_conda=args.no_use_conda,
                snakefile="run_anew.smk",
                snakemake_args=args.snakemake_args)

        run(cores=args.cores,
            no_conda=args.no_use_conda,
            snakemake_args=args.snakemake_args)
    except subprocess.CalledProcessError as e:
        print(f'Error in snakemake invocation: {e}', file=sys.stderr)
        sys.exit(e.returncode)
    sys.exit(0)


def run(
    cores: int = 4,
    no_conda: bool = False,
    snakefile: str = "Snakefile",
    workdir=None,
    snakemake_args: List[str] = None,
):
    with resource_path('blr', snakefile) as snakefile_path:
        cmd = ["snakemake", "-s", str(snakefile_path), "--cores", str(cores)]

        # Set defaults
        cmd += ["--printshellcmds"]

        if not no_conda:
            cmd += ["--use-conda"]

        if workdir is not None:
            cmd += ["--directory", str(workdir)]

        if "CONDA_ENVS" in os.environ:
            cmd += ["--conda-prefix", os.environ["CONDA_ENVS"]]

        if snakemake_args is not None:
            cmd += snakemake_args

        logger.debug(f"Command: {' '.join(cmd)}")
        subprocess.check_call(cmd)

"""
Run the BLR pipeline.

This is a small wrapper around Snakemake that sets some default parameters
"""
from importlib_resources import path as resource_path
import logging
import sys
import os

from snakemake import snakemake
from snakemake.utils import available_cpu_count

logger = logging.getLogger(__name__)


class SnakemakeError(Exception):
    pass


def add_arguments(parser):
    arg = parser.add_argument
    arg('--dryrun', '-n', default=False, action='store_true',
        help='Do not execute anything.')
    arg('--cores', '--jobs', '-j', metavar='N', type=int, default=available_cpu_count(),
        help='Run on at most N CPU cores in parallel. '
        'Default: %(default)s (all available cores).')
    arg('--keepgoing', '-k', default=False, action='store_true',
        help='If one job fails, finish the others.')
    arg('--unlock', default=False, action='store_true',
        help='Remove a lock on the working directory.')
    arg('--delete-all-output', default=False, action='store_true',
        help="Remove all files generated by the snakemake workflow. Use together with -n/--dry-run to list files "
             "without actually deleting anything. Write-protected files are not removed. Nevertheless, use with care! "
             "Default: %(default)s.")
    arg('--force-run', '-R', nargs='*',
        help="Force the re-execution or creation of the given rules or files. Use this option if you changed a rule "
             "and want to have all its output in your workflow updated. Default: %(default)s.")
    dags = parser.add_mutually_exclusive_group()
    dags.add_argument(
        "--dag", default=False, action='store_true',
        help="Print the dag in the graphviz dot language (requires graphviz to be installed). Default: %(default)s. "
             "To get output to pdf file, pipe output into dot as follows: blr run --dag | dot -Tpdf > dag.pdf."
    )
    dags.add_argument(
        "--filegraph", default=False, action='store_true',
        help="Print the file graph showing input/output file from rules in the graphviz dot language (requires "
             "graphviz to be installed). Default: %(default)s. To get output to pdf file, pipe output into dot "
             "as follows: blr run --filegraph | dot -Tpdf > filegraph.pdf."
    )
    arg("--snakemake-kws", nargs=2, metavar=("KEY", "VALUE"), action="append", default=[],
        help=f"Additional snakemake arguments not yet added to {__name__}. See the Snakemake API ("
             "https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html#) for options. Note that some "
             "options may still not be available.")
    arg('targets', nargs='*', default=[],
        help="File(s) to create. If omitted, the full pipeline is run. Include 'anew' if initializing from a "
             "previous analysis run(s).")


def main(args):

    try:
        if "anew" in args.targets:
            run(dryrun=args.dryrun,
                cores=args.cores,
                keepgoing=args.keepgoing,
                unlock=args.unlock,
                delete_all_output=args.delete_all_output,
                force_run=args.force_run,
                printdag=args.dag,
                printfilegraph=args.filegraph,
                snake_kws=dict(args.snakemake_kws),
                targets=None,
                snakefile="run_anew.smk")
            args.targets.remove("anew")

        run(dryrun=args.dryrun,
            cores=args.cores,
            keepgoing=args.keepgoing,
            unlock=args.unlock,
            delete_all_output=args.delete_all_output,
            force_run=args.force_run,
            printdag=args.dag,
            printfilegraph=args.filegraph,
            snake_kws=dict(args.snakemake_kws),
            targets=args.targets)
    except SnakemakeError:
        sys.exit(1)
    sys.exit(0)


def run(
    dryrun: bool = False,
    cores: int = 4,
    keepgoing: bool = False,
    unlock: bool = False,
    delete_all_output: bool = False,
    force_run=None,
    printdag: bool = False,
    printfilegraph: bool = False,
    snake_kws=None,
    targets=None,
    workdir=None,
    snakefile="Snakefile"
):
    snake_kws = {} if snake_kws is None else snake_kws

    conda_prefix = get_conda_prefix(snake_kws)

    # snakemake sets up its own logging, and this cannot be easily changed
    # (setting keep_logger=True crashes), so remove our own log handler
    # for now
    logger.root.handlers = []
    with resource_path('blr', snakefile) as snakefile_path:
        success = snakemake(
            snakefile_path,
            snakemakepath='snakemake',
            dryrun=dryrun,
            cores=cores,
            keepgoing=keepgoing,
            unlock=unlock,
            delete_all_output=delete_all_output,
            forcerun=force_run,
            printshellcmds=True,
            printdag=printdag,
            printfilegraph=printfilegraph,
            targets=targets,
            workdir=workdir,
            use_conda=True,
            printreason=dryrun,
            conda_prefix=conda_prefix,
            log_handler=[print_log_on_error],
            **snake_kws
        )
    if not success:
        raise SnakemakeError()


def print_log_on_error(msg):
    """Prints logs of failed rules in case of error"""
    if msg["level"] == "job_error" and msg["log"]:
        for log in msg["log"]:
            head = f"=== Output from log: '{log}' ==="
            print(head)
            if log.exists:
                with open(log) as f:
                    print(f.read().strip())
            print("-"*len(head))


def get_conda_prefix(snakemake_kws):
    if "conda_prefix" in snakemake_kws:
        return snakemake_kws["conda_prefix"]
    elif "CONDA_ENVS" in os.environ:
        return os.environ["CONDA_ENVS"]

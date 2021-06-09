"""
BLR is a pipeline for processing barcoded long reads
"""
import sys
import logging
import pkgutil
import importlib
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import blr.cli as cli_package
from blr import __version__

logger = logging.getLogger(__name__)


def main(commandline_arguments=None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(module)s - %(levelname)s: %(message)s")
    parser = ArgumentParser(description=__doc__, prog="blr")
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("--debug", action="store_true", default=False, help="Print debug messages")
    parser.add_argument("--profile", action="store_true", default=False,
                        help="Save profile info to blr_<subcommand>.prof")
    subparsers = parser.add_subparsers()

    # Import each module that implements a subcommand and add a subparser for it.
    # Each subcommand is implemented as a module in the cli subpackage.
    # It needs to implement an add_arguments() and a main() function.
    modules = pkgutil.iter_modules(cli_package.__path__)
    for _, module_name, _ in modules:
        module = importlib.import_module("." + module_name, cli_package.__name__)
        help = module.__doc__.strip().split("\n", maxsplit=1)[0]
        subparser = subparsers.add_parser(module_name, help=help, description=module.__doc__,
                                          formatter_class=RawDescriptionHelpFormatter)
        subparser.set_defaults(module=module)
        module.add_arguments(subparser)

    args = parser.parse_args(commandline_arguments)

    if args.debug:
        root = logging.getLogger()
        root.setLevel(logging.DEBUG)

    if not hasattr(args, "module"):
        parser.error("Please provide the name of a subcommand to run")
    else:
        module = args.module
        subcommand = module.main
        del args.module
        del args.debug
        profile = args.profile
        del args.profile

        # Print settings for module
        module_name = module.__name__.split('.')[-1]
        sys.stderr.write(f"SETTINGS FOR: {module_name} (version: {__version__})\n")
        for object_variable, value in vars(args).items():
            sys.stderr.write(f" {object_variable}: {value}\n")

        if profile:
            import cProfile
            profile_file = f'blr_{module_name}.prof'
            cProfile.runctx("subcommand(args)", globals(), dict(subcommand=subcommand, args=args),
                            filename=profile_file)
            logger.info(f"Writing profiling stats to '{profile_file}'.")
        else:
            subcommand(args)

    return 0


if __name__ == "__main__":
    sys.exit(main())

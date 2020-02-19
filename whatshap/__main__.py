import sys
import pkgutil
import importlib
import logging

import whatshap.cli as cli_package
from . import __version__
from .args import HelpfulArgumentParser
from .cli import CommandLineError


logger = logging.getLogger(__name__)


class NiceFormatter(logging.Formatter):
    """
    Do not prefix "INFO:" to info-level log messages (but do it for all other
    levels).

    Based on http://stackoverflow.com/a/9218261/715090 .
    """

    def format(self, record):
        if record.levelno != logging.INFO:
            record.msg = "{}: {}".format(record.levelname, record.msg)
        return super().format(record)


def setup_logging(debug):
    """
    Set up logging. If debug is True, then DEBUG level messages are printed.
    """
    handler = logging.StreamHandler()
    handler.setFormatter(NiceFormatter())
    root = logging.getLogger()
    root.addHandler(handler)
    root.setLevel(logging.DEBUG if debug else logging.INFO)


def ensure_pysam_version():
    from pysam import __version__ as pysam_version
    from distutils.version import LooseVersion

    if LooseVersion(pysam_version) < LooseVersion("0.8.1"):
        sys.exit("WhatsHap requires pysam >= 0.8.1")


def main(argv=sys.argv[1:]):
    ensure_pysam_version()
    parser = HelpfulArgumentParser(description=__doc__, prog="whatshap")
    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)
    parser.add_argument("--debug", action="store_true", default=False, help="Print debug messages")
    subparsers = parser.add_subparsers()

    # Import each module that implements a subcommand and add a subparser for it.
    # Each subcommand is implemented as a module in the cli subpackage.
    # It needs to implement an add_arguments() and a main() function.
    modules = pkgutil.iter_modules(cli_package.__path__)
    for _, module_name, _ in modules:
        module = importlib.import_module("." + module_name, cli_package.__name__)
        subparser = subparsers.add_parser(
            module_name,
            help=module.__doc__.strip().split("\n", maxsplit=1)[0],
            description=module.__doc__,
        )
        subparser.set_defaults(module=module, subparser=subparser)
        module.add_arguments(subparser)

    args = parser.parse_args(argv)
    setup_logging(args.debug)

    if not hasattr(args, "module"):
        parser.error("Please provide the name of a subcommand to run")
    else:
        module = args.module
        if hasattr(args.module, "validate"):
            subparser = args.subparser
            args.module.validate(args, subparser)
        del args.subparser
        del args.module
        del args.debug
        try:
            module.main(args)
        except CommandLineError as e:
            logger.error("whatshap error: %s", str(e))
            sys.exit(1)


if __name__ == "__main__":
    main()

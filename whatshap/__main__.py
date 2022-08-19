import ast
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
            record.msg = f"{record.levelname}: {record.msg}"
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


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    subcommand_name = get_subcommand_name(argv)
    module = importlib.import_module("." + subcommand_name, cli_package.__name__)

    parser = HelpfulArgumentParser(description=__doc__, prog="whatshap")
    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)
    parser.add_argument("--debug", action="store_true", default=False, help="Print debug messages")
    subparsers = parser.add_subparsers()
    subparser = subparsers.add_parser(
        subcommand_name,
        help=module.__doc__.strip().split("\n", maxsplit=1)[0],
        description=module.__doc__,
    )
    module.add_arguments(subparser)
    args = parser.parse_args(argv)
    setup_logging(args.debug)

    if hasattr(module, "validate"):
        module.validate(args, subparser)
    del args.debug
    try:
        module.main(args)
    except CommandLineError as e:
        logger.error("whatshap error: %s", str(e))
        logger.debug("Command line error. Traceback:", exc_info=True)
        sys.exit(1)


def get_subcommand_name(arguments) -> str:
    """
    Parse arguments to find out which subcommand was requested.

    This sets up a minimal ArgumentParser with the correct help strings.

    Because help is obtained from a module’s docstring, but importing each module
    makes startup slow, the modules are only parsed with the ast module and
    not fully imported at this stage.

    Return:
        subcommand name
    """
    parser = HelpfulArgumentParser(description=__doc__, prog="whatshap")
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers()

    for module_name, docstring in cli_modules(cli_package):
        help = docstring.strip().split("\n", maxsplit=1)[0].replace("%", "%%")
        subparser = subparsers.add_parser(
            module_name, help=help, description=docstring, add_help=False
        )
        subparser.set_defaults(module_name=module_name)
    args, _ = parser.parse_known_args(arguments)
    module_name = getattr(args, "module_name", None)
    if module_name is None:
        parser.error("Please provide the name of a subcommand to run")
    return module_name


def cli_modules(package):
    """
    Yield (module_name, docstring) tuples for all modules in the given package.
    """
    modules = pkgutil.iter_modules(package.__path__)
    for module in modules:
        spec = importlib.util.find_spec(package.__name__ + "." + module.name)
        with open(spec.origin) as f:
            mod_ast = ast.parse(f.read())
        docstring = ast.get_docstring(mod_ast, clean=False)
        yield module.name, docstring


if __name__ == "__main__":
    main()

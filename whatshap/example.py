"""
Print the square of a number

This is an example subcommand. The first line above is a short summary description
of what the command does. It is shown when running whatshap --help. Make it short!

To add a new subcommand:

* Copy this module to a new file and adjust it to your needs. The name of the
  module is identical to the subcommand name.
* Add the name of this module to the ``COMMANDS`` list in ``whatshap/__main__.py``.

The module must have add_arguments() and main() methods. The validate() method is optional and
can be used to validate the command-line options (use parser.error to raise an
error).
"""
import logging

logger = logging.getLogger(__name__)


def add_arguments(parser):
	add = parser.add_argument
	add('number', type=int, help='The number to square (at least 10)')


def validate(args, parser):
	if args.number <= 10:
		parser.error("Sorry, this is too simple!")


def main(args):
	print('{} * {} = {}'.format(args.number, args.number, args.number * args.number))

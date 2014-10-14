"""
Following Cython's recommendations: If pre-generated .c extension files are
found, then Cython is not run, even if it is installed.
Pass --cython on the command line to force Cython to be run.

If .c files are not found, such as in a fresh Git checkout, then Cython is
always run.
"""
import sys
import os.path
from distutils.core import setup, Extension
from distutils.version import LooseVersion

from whatshap import __version__

MIN_CYTHON_VERSION = '0.15'

if sys.version_info < (3, 3):
	sys.stdout.write("At least Python 3.3 is required.\n")
	sys.exit(1)


use_cython = False  # TODO not os.path.exists('whatshap/_somefile.c')

if '--cython' in sys.argv:
	use_cython = True
	sys.argv.remove('--cython')

# Try to find out whether a recent enough Cython is installed.
# If it is not, fall back to using the pre-compiled C sources.
# Pre-compiled sources are available only in the official releases,
# not within the Git repository.
if use_cython:
	try:
		from Cython import __version__ as cyversion
	except ImportError:
		sys.stdout.write(
			"ERROR: Cython is not installed. Install at least Cython >= " + str(MIN_CYTHON_VERSION) +
		    " to continue.\n")
		sys.exit(1)
	if LooseVersion(cyversion) < LooseVersion(MIN_CYTHON_VERSION):
		sys.stdout.write(
			"Error: Your Cython is at version '" + str(cyversion) +
			"', but at least version " + str(MIN_CYTHON_VERSION) + " is required.\n")
		sys.exit(1)

	from Cython.Build import cythonize

ext = '.pyx' if use_cython else '.c'

# TODO add your extensions here
extensions = [
	# Extension('whatshap._myextension', sources=['whatshap/_myextension' + ext]),
]
if use_cython:
	extensions = cythonize(extensions)

setup(
	name = 'whatshap',
	version = __version__,
	author = '',
	author_email = '',
	url = '',
	description = '',
	license = '',
	ext_modules = extensions,
	packages = ['whatshap'],
	#scripts = ['bin/...'],
	classifiers = [
		"Development Status :: 2 - Pre-Alpha",
		"Environment :: Console",
		"Intended Audience :: Science/Research",
		#"License :: OSI Approved :: MIT License", ???
		"Natural Language :: English",
		#"Programming Language :: Cython",
		"Programming Language :: Python",
		"Programming Language :: Python :: 3",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	]
)

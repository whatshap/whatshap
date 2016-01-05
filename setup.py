"""
While we use Cython as programming language for the extension modules, we
follow Cython's recommendation to distribute pre-generated .c/.cpp files.
Thus, Cython does not need to be installed on the machine of the user installing
WhatsHap. Pass --cython on the command line to force Cython to be run.

If the .c/.cpp files are not found or out of date, such as in a fresh Git
checkout, then Cython is always run.
"""
import sys
import os.path
from setuptools import setup, Extension
from distutils.version import LooseVersion

# set __version__
exec(next(open('whatshap/__init__.py')))

MIN_CYTHON_VERSION = '0.17'

if sys.version_info < (3, 2):
	sys.stdout.write("At least Python 3.2 is required.\n")
	sys.exit(1)
if sys.version_info < (3, 3):
	# need backport of contextlib.ExitStack
	extra_install_requires = ['contextlib2']
else:
	extra_install_requires = []


def out_of_date(extensions):
	"""
	Check whether any pyx source is newer than the corresponding generated
	C(++) source or whether any C(++) source is missing.
	"""
	for extension in extensions:
		for pyx in extension.sources:
			path, ext = os.path.splitext(pyx)
			if ext not in ('.pyx', '.py'):
				continue
			csource = path + ('.cpp' if extension.language == 'c++' else '.c')
			# When comparing modification times, allow five seconds slack:
			# If the installation is being run from pip, modification
			# times are not preserved and therefore depends on the order in
			# which files were unpacked.
			if not os.path.exists(csource) or (
				os.path.getmtime(pyx) > os.path.getmtime(csource) + 5):
				return True
	return False


def no_cythonize(extensions, **_ignore):
	"""
	Change file extensions from .pyx to .c or .cpp.

	Copied from Cython documentation
	"""
	for extension in extensions:
		sources = []
		for sfile in extension.sources:
			path, ext = os.path.splitext(sfile)
			if ext in ('.pyx', '.py'):
				if extension.language == 'c++':
					ext = '.cpp'
				else:
					ext = '.c'
				sfile = path + ext
			sources.append(sfile)
		extension.sources[:] = sources
	return extensions


def cythonize_if_necessary(extensions):
	if '--cython' in sys.argv:
		sys.argv.remove('--cython')
	elif out_of_date(extensions):
		sys.stdout.write('At least one C source file is missing or out of date.\n')
	else:
		return no_cythonize(extensions)

	try:
		from Cython import __version__ as cyversion
	except ImportError:
		sys.stdout.write(
			"ERROR: Cython is not installed. Install at least Cython version " +
			str(MIN_CYTHON_VERSION) + " to continue.\n")
		sys.exit(1)
	if LooseVersion(cyversion) < LooseVersion(MIN_CYTHON_VERSION):
		sys.stdout.write(
			"ERROR: Your Cython is at version '" + str(cyversion) +
			"', but at least version " + str(MIN_CYTHON_VERSION) + " is required.\n")
		sys.exit(1)

	from Cython.Build import cythonize
	return cythonize(extensions)

extensions = [
	Extension('whatshap._core',
		sources=['whatshap/_core.pyx',
			'src/columncostcomputer.cpp', 'src/columnindexingiterator.cpp',
			'src/columnindexingscheme.cpp', 'src/dptable.cpp',
			'src/entry.cpp', 'src/graycodes.cpp', 'src/read.cpp',
			'src/readset.cpp', 'src/columniterator.cpp', 'src/indexset.cpp'
		], language='c++', extra_compile_args=["-std=c++0x"],),
]
extensions = cythonize_if_necessary(extensions)

setup(
	name = 'whatshap',
	version = __version__,
	author = 'WhatsHap authors',
	author_email = 'whatshap@cwi.nl',
	url = 'https://bitbucket.org/whatshap/whatshap/',
	description = 'phase genomic variants using DNA sequencing reads',
	license = 'MIT',
	ext_modules = extensions,
	packages = ['whatshap', 'whatshap.scripts'],
	scripts = ['bin/whatshap', 'bin/phasingstats'],
	install_requires = ['pysam>=0.8.1', 'PyVCF'] + extra_install_requires,
	classifiers = [
		"Development Status :: 4 - Beta",
		"Environment :: Console",
		"Intended Audience :: Science/Research",
		"License :: OSI Approved :: MIT License",
		"Natural Language :: English",
		"Programming Language :: Cython",
		"Programming Language :: Python",
		"Programming Language :: Python :: 3",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	]
)

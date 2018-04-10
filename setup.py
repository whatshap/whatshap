"""
While we use Cython as programming language for the extension modules, we
follow Cython's recommendation to distribute pre-generated .c/.cpp files.
Thus, Cython does not need to be installed on the machine of the user installing
WhatsHap.
"""
import sys
import os.path
from setuptools import setup, Extension
from distutils.version import LooseVersion
from distutils.command.sdist import sdist as _sdist
from distutils.command.build_ext import build_ext as _build_ext
from distutils.sysconfig import customize_compiler
import versioneer


MIN_CYTHON_VERSION = '0.17'

if sys.version_info < (3, 4):
	sys.stdout.write("At least Python 3.4 is required.\n")
	sys.exit(1)


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


def check_cython_version():
	"""exit if Cython not found or out of date"""
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


def CppExtension(name, sources):
	return Extension(name, sources=sources, language='c++',
		extra_compile_args=['-std=c++11'], undef_macros = [ "NDEBUG" ])


extensions = [
	CppExtension('whatshap.core',
		sources=['whatshap/core.pyx',
			'src/pedigree.cpp',
			'src/pedigreedptable.cpp', 'src/pedigreecolumncostcomputer.cpp',
			'src/columnindexingiterator.cpp', 'src/columnindexingscheme.cpp',
			'src/entry.cpp', 'src/graycodes.cpp', 'src/read.cpp',
			'src/readset.cpp', 'src/columniterator.cpp', 'src/indexset.cpp',
			'src/pedigreepartitions.cpp', 'src/phredgenotypelikelihoods.cpp',
			'src/genotyper.cpp', 'src/genotypedistribution.cpp',
			'src/genotypedptable.cpp', 'src/genotypecolumncostcomputer.cpp',
			'src/backwardcolumniterator.cpp', 'src/transitionprobabilitycomputer.cpp'
		]),
	CppExtension('whatshap.priorityqueue', sources=['whatshap/priorityqueue.pyx']),
	CppExtension('whatshap.align', sources=['whatshap/align.pyx']),
	CppExtension('whatshap._variants', sources=['whatshap/_variants.pyx']),
]

cmdclass = versioneer.get_cmdclass()


class build_ext(cmdclass.get('build_ext', _build_ext)):
	def run(self):
		# If we encounter a PKG-INFO file, then this is likely a .tar.gz/.zip
		# file retrieved from PyPI that already includes the pre-cythonized
		# extension modules, and then we do not need to run cythonize().
		if os.path.exists('PKG-INFO'):
			no_cythonize(extensions)
		else:
			# Otherwise, this is a 'developer copy' of the code, and then the
			# only sensible thing is to require Cython to be installed.
			check_cython_version()
			from Cython.Build import cythonize
			self.extensions = cythonize(self.extensions)
		super().run()

	def build_extensions(self):
		# Remove the warning about “-Wstrict-prototypes” not being valid for C++,
		# see http://stackoverflow.com/a/36293331/715090
		customize_compiler(self.compiler)
		if self.compiler.compiler_so[0].endswith('clang'):
			# Clang needs this option in order to find the unordered_set header
			print('detected clang, using option -stdlib=libc++')
			self.compiler.compiler_so.append('-stdlib=libc++')
		try:
			self.compiler.compiler_so.remove("-Wstrict-prototypes")
		except (AttributeError, ValueError):
			pass
		super().build_extensions()


class sdist(cmdclass.get('sdist', _sdist)):
	def run(self):
		# Make sure the compiled Cython files in the distribution are up-to-date
		from Cython.Build import cythonize
		check_cython_version()
		cythonize(extensions)
		super().run()


cmdclass['build_ext'] = build_ext
cmdclass['sdist'] = sdist

with open('doc/README.rst', encoding='utf-8') as f:
	long_description = f.read()

if sys.version_info < (3, 5):
	requires = ['typing']
else:
	requires = []

setup(
	name = 'whatshap',
	version = versioneer.get_version(),
	author = 'WhatsHap authors',
	author_email = 'whatshap@cwi.nl',
	url = 'https://bitbucket.org/whatshap/whatshap/',
	description = 'phase genomic variants using DNA sequencing reads',
	long_description = long_description,
	license = 'MIT',
	cmdclass = cmdclass,
	ext_modules = extensions,
	packages = ['whatshap'],
	entry_points={'console_scripts': ['whatshap = whatshap.__main__:main']},
	install_requires = ['pysam>=0.14.0', 'PyVCF', 'pyfaidx', 'xopen'] + requires,
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

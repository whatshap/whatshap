"""Misc. tests"""
import whatshap.__main__


def test_main():
	"""Test the main command-line module"""
	try:
		whatshap.__main__.main(['--version'])
	except SystemExit as e:
		if e.code != 0:
			raise


def test_pysam_version():
	from pysam import __version__ as pysam_version
	from distutils.version import LooseVersion
	assert LooseVersion(pysam_version) >= LooseVersion("0.8.1")

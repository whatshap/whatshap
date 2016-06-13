"""Test the main command-line module"""
import whatshap.__main__


def testmain():
	try:
		whatshap.__main__.main(['--version'])
	except SystemExit as e:
		if e.code != 0:
			raise

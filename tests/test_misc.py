"""Misc. tests"""
import whatshap.__main__


def test_main():
    """Test the main command-line module"""
    try:
        whatshap.__main__.main(["--version"])
    except SystemExit as e:
        if e.code != 0:
            raise

__all__ = ["__version__"]

from ._version import version as __version__


def _load_cpp_classes():
    import sys
    import os

    flags = sys.getdlopenflags()
    sys.setdlopenflags(os.RTLD_GLOBAL | os.RTLD_NOW)
    import whatshap.core  # noqa
    import whatshap.polyphase.solver  # noqa

    sys.setdlopenflags(flags)


_load_cpp_classes()

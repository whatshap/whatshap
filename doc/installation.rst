.. _installation:

============
Installation
============

Installation of WhatsHap is easiest if you use Conda.


Installation with Conda
-----------------------

First, ensure you have Conda (miniconda or Anaconda) installed and made the
proper settings to enable the “bioconda” channel. For that, follow
`the bioconda installation instructions <https://bioconda.github.io/#install-conda>`_.

Then WhatsHap can be installed with this command::

    conda install whatshap

If you have Conda, but not enabled bioconda, use this command::

    conda install -c bioconda -c conda-forge whatshap


Installation with pip
---------------------

Before you can `pip install`, you need to install dependencies that pip cannot
install for you. WhatsHap is implemented in C++ and Python. You need to have a
C++ compiler, Python 3.3 (or later) and the corresponding Python header files.
Python 3.5 is slightly faster than 3.4. In Ubuntu, installing the packages
``build-essential`` and ``python3-dev`` will take care of all required
dependencies.

WhatsHap can then be installed with pip::

	pip3 install --user whatshap

This installs WhatsHap into ``$HOME/.local/bin``.  Then add
``$HOME/.local/bin`` to your ``$PATH`` and run the tool::

    export PATH=$HOME/.local/bin:$PATH
    whatshap --help

Alternatively, you can also install WhatsHap into a virtual environment if you
are familiar with that.


Installing an unreleased development version
--------------------------------------------

If you for some reason want to use the most recent development version of
WhatsHap, you can install it in the following way. These instructions will
create a virtual environment in the directory ``whatshap-env`` that contains
WhatsHap. Simply delete that directory to uninstall the software. Other WhatsHap
versions you may have installed in other locations remain unaffected. ::

	python3 -m venv whatshap-env
	whatshap-env/bin/pip install Cython
	whatshap-env/bin/pip install git+https://bitbucket.org/whatshap/whatshap

You can then run WhatsHap like this::

	whatshap-env/bin/whatshap --version

You should see a version number like ``0.17+103.g71e5b3c``, which means that
this version is 103 Git commits ahead of version 0.17.

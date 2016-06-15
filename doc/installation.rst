============
Installation
============


Requirements
------------

WhatsHap is implemented in C++ and Python. You need to have a C++ compiler,
Python 3.3 (or later) and the corresponding Python header files. Python 3.5
is slightly faster than 3.4. In Ubuntu, make sure the packages
``build-essential`` and ``python3-dev`` are installed.


Installation
------------

WhatsHap can be installed with pip::

	pip3 install --user whatshap

This installs WhatsHap into ``$HOME/.local/bin``.  Then add
``$HOME/.local/bin`` to your ``$PATH`` and run the tool::

	export PATH=$HOME/.local/bin:$PATH
	whatshap --help


Installation with Conda
-----------------------

WhatsHap is also available as a conda package from the `bioconda
channel <https://bioconda.github.io/>`_::

    conda install -c bioconda whatshap

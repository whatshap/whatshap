============
Installation
============


Requirements
------------

WhatsHap is implemented in C++ and Python. You need to have a C++ compiler,
Python 3.2 (or later) and the corresponding Python header files. In Ubuntu,
make sure the packages ``build-essential`` and ``python3-dev`` are installed.


Quickstart
----------

As soon as there is a release, this should work::

	pip3 install --user WhatsHap

Then add ``$HOME/.local/bin`` to your ``$PATH`` and run the tool::

	export PATH=$HOME/.local/bin:$PATH
	whatshap --help


Regular installation
--------------------

There is currently no release of WhatsHap, so you need to install it from the
Git repository instead. Make sure you also have installed Cython::

	pip3 install --user Cython
	pip3 install --user https://bitbucket.org/whatshap/whatshap/get/master.tar.gz

This installs WhatsHap into ``$HOME/.local/bin``. The Cython requirement will
be dropped when there is a first release.

You can also use a virtualenv instead, but you need to make sure that you have
installed Cython into the virtualenv before installing WhatsHap::

	virtualenv -p python3 venv  # Creates a virtualenv in the directory 'venv'
	venv/bin/pip3 install Cython   # Installs Cython into the virtualenv
	venv/bin/pip3 install https://bitbucket.org/whatshap/whatshap/get/master.tar.gz

If you get errors while installing Cython, try to add
``--install-option="--no-cython-compile"`` to the command, see also
`issue 43 <https://bitbucket.org/whatshap/whatshap/issue/43/>`_.


Development installation
------------------------

For development, make sure that you install Cython. We also recommend using a
virtualenv. This sequence of commands should work::

	git clone https://bitbucket.org/whatshap/whatshap
	cd whatshap
	virtualenv -p python3 venv
	venv/bin/pip3 install Cython nose
	venv/bin/python3 setup.py develop

Then you can run WhatsHap like this::

	venv/bin/whatshap --help

The tests can be run like this::

	venv/bin/nosetests

If you use a nosetests binary from your system, it is usually called
`nosetests3`.


Development installation (alternative)
--------------------------------------

Alternatively, if you do not want to use virtualenv, you can do the following::

	git clone https://bitbucket.org/whatshap/whatshap.git
	cd whatshap
	python3 setup.py build_ext -i --cython
	bin/whatshap

This requires Cython, pysam, and pyvcf to be installed.


Installing other Python versions in Ubuntu
------------------------------------------

Ubuntu comes with one default Python 3 version, and in order to test WhatsHap
with other Python versions (3.2, 3.3 and 3.4), use the “deadsnakes” repository.
Ensure you have the following packages::

	sudo apt-get install build-essential python-software-properties

Then get and install the desired Python versions. For example, for Python 3.2::

	sudo add-apt-repository ppa:fkrull/deadsnakes
	sudo apt-get update
	sudo apt-get install python3.2-dev python3-setuptools

If pip and virtualenv are not available, install them (Since they are so essential,
we use sudo to install them system-wide, but you can also install them into
your $HOME by omitting the sudo and adding the ``--user`` option instead)::

	sudo easy_install3 pip
	sudo pip3 install virtualenv

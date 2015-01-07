============
Installation
============


Requirements
------------

WhatsHap is implemented in C++ and Python. You need to have a C++ compiler,
Python 3.3 (or later) and the corresponding Python header files. In Ubuntu,
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
Bitbucket repository instead. Make sure you also have installed Cython::

	pip3 install --user Cython
	pip3 install --user https://bitbucket.org/whatshap/whatshap/get/master.tar.gz

This installs WhatsHap into ``$HOME/.local/bin``. The Cython requirement will
be dropped when there is a first release.

You can also use a virtualenv instead, but you need to make sure that you have
installed Cython into the virtualenv before installing WhatsHap::

	virtualenv -p python3 venv
	venv/bin/pip3 install Cython
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
	venv/bin/pip3 install Cython
	venv/bin/python3 setup.py develop

Then you can run WhatsHap like this::

	venv/bin/whatshap --help

Development installation (alternative)
--------------------------------------
Alternatively, if you do not want to use virtualenv, you can do the following::

	git clone https://bitbucket.org/whatshap/whatshap.git
	cd whatshap
	python3 setup.py build_ext -i --cython
	./bin/whatshap

This requires Cython, pysam, and pyvcf to be installed.

Ubuntu 12.04 LTS
----------------

This Ubuntu release does not have a recent enough Python, so you need to install
it first. From a freshly installed Ubuntu 12.04, you would need to do this
first::

	sudo apt-get install build-essential python-software-properties

Then get and install Python 3.3 (or choose 3.4)::

	sudo add-apt-repository ppa:fkrull/deadsnakes
	sudo apt-get update
	sudo apt-get install python3.3-dev python3-setuptools

Then install pip and virtualenv. (Since they are so essential, we use sudo to
install them system-wide, but you can also install them into your $HOME by
omitting the sudo and adding the ``--user`` option instead)::

	sudo easy_install3 pip
	sudo pip3 install virtualenv

Then continue with the regular or development instructions above to install
WhatsHap.

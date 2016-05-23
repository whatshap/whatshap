Developing
==========

This is the developer documentation for WhatsHap.


Development installation
------------------------

For development, make sure that you install Cython and tox. We also recommend
using a virtualenv. This sequence of commands should work::

	git clone https://bitbucket.org/whatshap/whatshap
	cd whatshap
	virtualenv -p python3 venv
	venv/bin/pip3 install Cython nose tox
	venv/bin/pip3 install -e .

Then you can run WhatsHap like this (or activate the virtualenv and omit the
``venv/bin`` part)::

	venv/bin/whatshap --help

The tests can be run like this::

	venv/bin/tox


Development installation (alternative)
--------------------------------------

Alternatively, if you do not want to use virtualenv, you can do the following::

	git clone https://bitbucket.org/whatshap/whatshap.git
	cd whatshap
	python3 setup.py build_ext -i
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


Debugging
---------

	$ gdb python3
	(gdb) run -m nose

After you get a SIGSEGV, let gdb print a backtrace:

	(gdb) bt

============
Installation
============

There is currently no release of WhatsHap, so you need to install it from the
Bitbucket repository::

	pip3 install --user Cython
	pip3 install --user https://bitbucket.org/marcelm/whatshap/get/master.tar.gz

This installs WhatsHap into ``$HOME/.local/bin``.

You can also use a virtualenv, but you need to make sure that you have installed
Cython into the virtualenv before installing WhatsHap::

	virtualenv -p python3 venv
	venv/bin/pip install Cython
	pip3 install https://bitbucket.org/marcelm/whatshap/get/master.tar.gz

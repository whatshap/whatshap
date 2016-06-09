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

Here is one way to get a backtrace from gdb (assuming the bug occurs while
running the tests)::

	$ gdb python3
	(gdb) run -m nose

After you get a SIGSEGV, let gdb print a backtrace:

	(gdb) bt


Making a release
----------------

If this is the first time you attempt to upload a distribution to PyPI, create a
configuration file named ``.pypirc`` in your home directory with the following
contents::

	[distutils]
	index-servers =
	    pypi

	[pypi]
	username=my-user-name
	password=my-password

See also `this blog post about getting started with
PyPI <http://peterdowns.com/posts/first-time-with-pypi.html>`_. In particular,
note that a ``%`` in your password needs to be doubled and that the password
must *not* be put between quotation marks even if it contains spaces.

#. Ensure all tests pass.

#. Update the ``__version__`` in ``whatshap/__init__.py``.

#. Make sure the changelog is up-to-date. Set also the correct version number in the changelog.

#. Commit the version change.

#. Create a tag::

       git tag v0.1

#. Create and upload the distribution (``.tar.gz`` file) to PyPI::

       python3 setup.py sdist upload

#. Push, uploading the tag::

       git push --tags

   If something went wrong, increment the revision in the version number and re-do.
   That is, go from version x.y to x.y.1. Make sure that the tag in the repository
   always corresponds to the tarball that is on PyPI.

#. Update the `bioconda recipe <https://github.com/bioconda/bioconda-recipes/blob/master/recipes/whatshap/meta.yaml>`_.
   It is probly easiest to edit the recipe via the web interface and send in a
   pull request. Ensure that the list of dependencies (the ``requirements:``
   section in the recipe) is in sync with the ``setup.py`` file.

   Since this is just a version bump, the pull request does not need a
   review by other bioconda developers. As soon as the tests pass and if you
   have the proper permissions, it can be merged directly.


Adding a new subcommand
-----------------------

Follow the instructions in ``whatshap/example.py``.

Developing
==========

The `WhatsHap source code is on GitHub <https://github.com/whatshap/whatshap/>`_.
WhatsHap is developed in Python 3, Cython and C++.


Development installation
------------------------

For development, make sure that you install Cython and tox. We also recommend
using a virtualenv. This sequence of commands should work (use
``https://github.com/whatshap/whatshap.git`` as URL if you do not have a
GitHub account)::

    git clone git@github.com:whatshap/whatshap.git
    cd whatshap
    python3 -m venv venv
    source venv/bin/activate
    pip install -e .[dev]

The last command installs also all the development dependencies.
Omit the ``[dev]`` to leave them out.

Next, you can run WhatsHap like this::

    whatshap --help


Development installation when using Conda
-----------------------------------------

If you are familiar with `Conda <https://docs.conda.io/en/latest/>`_, you can
also use a Conda environment for developing WhatsHap. We recommend that you
use Conda only to install Python itself and let the rest of the dependencies
be handled by ``pip``::

    conda create -n whatshap-dev python=3.10
    conda activate whatshap-dev
    pip install -e .[dev]


Running tests
-------------

While in the virtual (or Conda) environment, you can run the tests for the
current Python version like this::

    pytest

Whenever you change any Cython code (``.pyx`` files), you need to re-run the
``pip install -e .`` step in order to compile it.

Optionally, to run tests for all supported Python versions, you can run
`tox <https://tox.readthedocs.io/>`_. It creates separate virtual environments for each Python
version, installs WhatsHap into it and then runs the tests. It also tests documentation generation
with ``sphinx``. Run it like this::

    tox

If ``tox`` is installed on the system, you do not need to be inside a virtual environment for this.
However, you need to have all tested Python versions installed on the system! See the instructions
below for how to do this on Ubuntu.


Code style
----------

Python code needs to be formatted with `Black <https://github.com/psf/black>`_.
Either run `black whatshap tests` manually before you commit or use the
`pre-commit <https://pre-commit.com/>`_ framework to automate this.


Installing other Python versions in Ubuntu
------------------------------------------

Ubuntu comes with one default Python 3 version, and in order to test WhatsHap
with older or newer Python versions, follow the instructions for enabling the
`“deadsnakes” repository <https://launchpad.net/~deadsnakes/+archive/ubuntu/ppa>`_.
After you have done so, ensure you have the following packages::

    sudo apt install build-essential python-software-properties

Then get and install the desired Python versions. Make sure you install the ``-dev`` package.
For example, for Python 3.10::

    sudo apt update
    sudo apt install python3.10-dev


Debugging
---------

Here is one way to get a backtrace from gdb (assuming the bug occurs while
running the tests)::

    $ gdb python3
    (gdb) run -m pytest

After you get a SIGSEGV, let gdb print a backtrace:

    (gdb) bt

Another way is to set ``PYTHONFAULTHANDLER=1``::

    PYTHONFAULTHANDLER=1 pytest -vxs tests/test_run_whatshap.py


Wrapping C++ classes
--------------------

The WhatsHap phasing algorithm is written in C++, as are many of the core
data structures such as the “Read” class. To make the C++ classes usable from
Python, we use Cython to wrap the classes. All these definitions are spread
across multiple files. To add new attributes or methods to an existing class
or to add a new class, changes need to be made in different places.

Let us look at the “Read” class. The following places in the code may need to
be changed if the Read class is changed or extended:

* ``src/read.cpp``: Implementation of the class (C++).
* ``src/read.h``: Header with the class declaration (also normal C++).
* ``whatshap/cpp.pxd``: Cython declarations of the class. This repeats – using
  the Cython syntax this time – a subset of the information from the
  ``src/read.h`` file. This duplication is required because Cython
  cannot read ``.h`` files (it would need a full C++ parser for that).

  Note that the ``cpp.pxd`` file contains definitions for *all* the ``.h``
  headers. (It would be cleaner to have them in separate ``.pxd`` files, but
  this leads to problems when linking the compiled files.)
* ``whatshap/core.pxd``: This contains declarations of all *Cython* classes
  wrapping C++ classes. Note that the class ``Read`` in this file has the
  same name as the C++ class, but that it is not the same as the C++ one!
  The distinction is made by prefixing the C++ class with ``cpp.``, which is
  the name of the module in which it is declared in (that is, the C++ class
  ``Read`` is declared in ``cpp.pxd``). The wrapping (Cython) class ``Read``
  stores the C++ class in an attribute named ``thisptr``. If you add a new
  class, it needs to be added to this file. If you only modify an existing one,
  you probably do not need to change this file.
* ``whatshap/core.pyx``: The Cython implementation of the wrapper classes.
  Again, the name ``Read`` by itself is the Python wrapper class and
  ``cpp.Read`` is the name for the C++ class.

Before adding yet more C++ code, which then requires extra code for wrapping it,
consider writing an implementation in Cython instead. See ``readselect.pyx``,
for example, which started out as a Python module and was then transferred to
Cython to make it faster. Here, the Cython code is not merely a wrapper, but
contains the implementation itself.


Writing documentation
---------------------

Documentation is located in the ``doc/`` subdirectory. It is written in
`reStructuredText format <http://docutils.sourceforge.net/docs/user/rst/quickref.html>`_
and is translated by `Sphinx <http://www.sphinx-doc.org/>`_ into HTML format.

Documentation is hosted on `Read the Docs <https://readthedocs.org/>`_.
It is built automatically whenever a commit is made. The documentation in the
``master`` branch should be visible at
`https://whatshap.readthedocs.io/en/latest/ <https://whatshap.readthedocs.io/en/latest/>`_
and documentation for the most recent released version should be visible at
`https://whatshap.readthedocs.io/en/stable/ <https://whatshap.readthedocs.io/en/stable/>`_.

To generate documentation locally, ensure that you installed sphinx and the
add-ons necessary to build documentation (running ``pip install -e .[docs]`` will
take care of this). Then go into the ``doc/`` directory and run ``make``. You can
then open ``doc/_build/html/index.html`` in your browser. The theme that is
used is a bit different from the one used on Read the Docs.


Making a release
----------------

#. Update ``CHANGES.rst``: Set the correct version number and ensure that
   all nontrivial, user-visible changes are listed.

#. Ensure you have no uncommitted changes in the working copy.

#. Run ``tox``.

#. Tag the current commit with the version number (there must be a ``v`` prefix)::

       git tag -a -m "Version 0.1" v0.1

   To release a development version, use a ``dev`` version number such as
   ``v0.17.dev1``.
   Users will only get these when they use ``pip install --pre``.

#. Push the tag::

       git push --tags

#. Wait for the GitHub Action to finish. It will deploy the sdist and wheels to
   PyPI if everything worked correctly.

#. Update the `Bioconda recipe <https://github.com/bioconda/bioconda-recipes/blob/master/recipes/whatshap/meta.yaml>`_. It is easiest to wait for the Bioconda bot to open a PR. Ensure
   that the list of dependencies (the ``requirements:``
   section in the recipe) is in sync with the ``setup.py`` file. If changes are
   necessary to the bot-generated PR, just add your own commits to that PR.

If something went wrong, fix the problem and follow the above instructions again,
but with an incremented revision in the version number. That is, go from version
x.y to x.y.1. PyPI will not allow you to change a version that has already been
uploaded.


Adding a new subcommand
-----------------------

Use one of the modules in ``whatshap/cli/`` as a template. All modules in
that directory are automatically used as subcommands.


Download count statistics
-------------------------

Some statistics for the PyPI package are available at
`pypistats.org <https://pypistats.org/packages/whatshap>`_.

Here is a query for Google BigQuery that shows download counts (from PyPI)
since a given date, broken down by version ::

    SELECT
        file.project,
        file.version,
        COUNT(*) as total_downloads,
    FROM
        TABLE_DATE_RANGE(
            [the-psf:pypi.downloads],
            TIMESTAMP("20170101"),
            CURRENT_TIMESTAMP()
        )
    WHERE
        file.project = 'whatshap'
    GROUP BY
        file.project, file.version

Statistics for the Conda package are available on the
`WhatsHap package detail page <https://anaconda.org/bioconda/whatshap/>`_.

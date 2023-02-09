# Sphinx configuration

import os
import sys
import time

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.insert(0, os.path.abspath(os.pardir))

extensions = [
    "sphinx.ext.mathjax",
    "sphinx_issues",
    "sphinx_better_subsection",
]

templates_path = ["_templates"]

source_suffix = ".rst"

master_doc = "index"

authors = "Murray Patterson, Alexander Schönhuth, Tobias Marschall, Marcel Martin and WhatsHap contributors"
project = "whatshap"
copyright = "2014-{}, {}".format(time.gmtime().tm_year, authors)

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
from pkg_resources import get_distribution

release = get_distribution("whatshap").version
# Read The Docs modifies the conf.py script and we therefore get
# version numbers like 0.12+0.g27d0d31
if os.environ.get("READTHEDOCS") == "True":
    version = ".".join(release.split(".")[:2])
else:
    version = release

release = version

issues_github_path = "whatshap/whatshap"

exclude_patterns = ["_build"]

# The reST default role (used for this markup: `text`) to use for all
# documents.
# default_role = None

pygments_style = "sphinx"

html_theme = "default"

html_static_path = ["_static"]

#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
import importlib

sys.path.insert(0, os.path.abspath("../.."))

project = "idaes-gtep"
copyright = "2018-2025, IDAES"
author = "IDAES Development Team"
release = importlib.metadata.version("gtep")
version = importlib.metadata.version("gtep")

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx_collections",
    "nbsphinx",
]

templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

on_rtd = os.environ.get("READTHEDOCS", None) == "True"
if not on_rtd:  # only import and set the theme if we're building docs locally
    import sphinx_rtd_theme

    # html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]


# sphinx_collections is used to include notebooks from the main gtep directory
# in the Sphinx documentation. When the docs are built, this will copy all
# notebooks into a _collections/notebooks directory where Sphinx can find them.
collections = {
    "notebooks": {
        "driver": "copy_folder",
        "source": "../../gtep/",
        "target": "notebooks/",
        "ignore": ["*.py"],
    }
}

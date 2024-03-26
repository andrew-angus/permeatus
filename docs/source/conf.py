# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------

project = 'Permeatus'
copyright = '2024, Andrew Angus'
author = 'Andrew Angus, Lukasz Figiel'

# The full version, including alpha/beta/rc tags
release = '1.0.0'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinxcontrib.bibtex',
    'autoapi.extension',
]

# Tell autoapi where package is
autoapi_dirs = ['../../permeatus']

# Tell autoapi to skip class attributes
def skippers(app, what, name, obj, skip, options):
    if what == "attribute":
       skip = True
    return skip

def setup(sphinx):
   sphinx.connect("autoapi-skip-member", skippers)

# Bibliography file
bibtex_bibfiles = ["Hydrogen.bib"]

# Add any paths that contain templates here, relative this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Custom css file
html_css_files = ['custom.css']

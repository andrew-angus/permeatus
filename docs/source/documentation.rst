Documentation
=============

.. role:: bash(code)
   :language: bash

Editing
-------

Documentation source files sit in :bash:`./docs/source`, with respect to the package root directory.
Inside is the Python configuration file :bash:`conf.py`, and several :bash:`.rst` files which,
represent separate pages of the documentation. The bibilography file :bash:`Hydrogen.bib` is also
contained here, and can be updated if new citations are desired.

Compiling
---------

The documentation requires the Sphinx package, and the bibtex and autoapi extensions, which
can be installed with:

.. code-block:: console

   $ pip install sphinxcontrib-bibtex
   $ pip install sphinx-autoapi

Compiling the documentation to HTML is then achieved by first ensuring that you are in the
:bash:`./docs` folder, and running:

.. code-block:: console

   $ make html

The resulting HTML files are located in :bash:`./docs/build/html`, and opening 
:bash:`./docs/build/html/index.html` with a browser will open the documentation home page.

Similarly, if Latex is installed on the device, then compiling to PDF is achieved with:

.. code-block:: console

   $ make latexpdf

The resulting PDF will be located in :bash:`./docs/build/latex/permeatus.pdf`.

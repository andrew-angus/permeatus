Usage
=====

Cloning Repository
------------------

To acquire the code, clone from GitHub:

.. code-block:: console

   $ git clone https://github.com/lwfig/permeatus.git

.. note::

   For now, the GitHub repository is private, and access must be requested by
   sending an email to: l.w.figiel@warwick.ac.uk


.. _installation:

Installation
------------

.. role:: bash(code)
   :language: bash

To use permeatus, first ensure you are in the root directory of the 
repository (:bash:`cd permeatus` after cloning), and install it using pip:

.. code-block:: console

   $ pip install .

Basic Use
---------

For a demonstration of how to use the software, a tutorial Jupyter notebook is available
at :bash:`./tutorial/tutorial.ipynb`, with respect to the package root directory.

Development
-----------

If it is desired to simulate a 2D system, only the mesh function must be added, in the
style of those in the homogenisation module. The rest of the codes functionality will
then be utilised automatically.

#!/bin/python3
# Author: Andrew Angus

from setuptools import setup,find_packages

with open("README.md", "r") as fh:
  long_description = fh.read()

setup(
  name="permeatus",
  version="0.0.1",

  author="Andrew Angus",
  author_email="andrew.c.angus@warwick.ac.uk",

  packages=find_packages(include=['permeatus','permeatus.*']),
  include_package_data=True,

  url="https://github.com/andrewanguswarwick/permeatus",

  description="Permeation modelling tools built around ABAQUS.",
  long_description=long_description,
  long_description_content_type="text/markdown",

  python_requires='>=3.9',
  install_requires=[
    "numpy",
    #"setuptools<66.0.0",
    "matplotlib",
    ],
)

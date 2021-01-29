#!/usr/bin/env python

from setuptools import setup
import glob

from ragtag_utilities.utilities import get_ragtag_version

with open("README.md", "r") as fh:
    long_description = fh.read()

scripts = glob.glob("*.p*")

version = get_ragtag_version()[1:]

setup(
    name='RagTag',
    version=version,
    author='Michael Alonge',
    author_email='malonge11@gmail.com',
    description='Fast reference-guided genome assembly scaffolding',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/malonge/RagTag",
    packages=['ragtag_utilities'],
    package_dir={'ragtag_utilities': 'ragtag_utilities/'},
    license="MIT",
    classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
    ],
    install_requires=[
              'intervaltree',
              'numpy',
              'pysam',
              'networkx',
          ],
    python_requires='>=3.6',
    scripts=scripts,
    zip_safe=True
)
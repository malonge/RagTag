#!/usr/bin/env python

from setuptools import setup
import glob

from ragtag_utilities.utilities import get_ragtag_version

scripts = glob.glob("*.p*")

version = get_ragtag_version()[1:]

setup(
    name='RagTag',
    version=version,
    description='A tool to order and orient genome assembly contigs via mappings to a reference genome.',
    author='Michael Alonge',
    author_email='malonge11@gmail.com',
    packages=['ragtag_utilities'],
    package_dir={'ragtag_utilities': 'ragtag_utilities/'},
    install_requires=[
              'intervaltree',
              'numpy',
              'pysam',
          ],
    scripts=scripts,
    zip_safe=True
)
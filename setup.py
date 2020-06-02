#!/usr/bin/env python

from setuptools import setup
import glob

scripts = glob.glob("*.p*")

setup(
    name='RagTag',
    version='1.0.0',
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
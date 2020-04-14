#!/usr/bin/env python

from setuptools import setup
import glob

scripts = glob.glob("*.p*")

setup(
    name='RaGOO2',
    version='2.0.0',
    description='A tool to order and orient genome assembly contigs via alignments to a reference genome.',
    author='Michael Alonge',
    author_email='malonge11@gmail.com',
    packages=['ragoo2_utilities'],
    package_dir={'ragoo2_utilities': 'ragoo2_utilities/'},
    install_requires=[
              'intervaltree',
              'numpy',
              'pysam',
          ],
    scripts=scripts,
    zip_safe=True
)
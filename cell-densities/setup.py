#!/usr/bin/env python

import imp
import sys

from setuptools import setup, find_packages

if sys.version_info < (3, 6):
    sys.exit("Sorry, Python < 3.6 is not supported")

VERSION = imp.load_source("", "nrrdhlp/version.py").__version__

setup(
    name="cell-densities",
    author=["Michael W. Reimann"],
    author_email="michael.reimann@epfl.ch",
    version=VERSION,
    description="validate cell densities",
    long_description="validate prescribed and built cell densities",
    url="http://bluebrain.epfl.ch",
    license="LGPL-3.0",
    install_requires=["numpy",
                      "pandas",
                      "voxcell",
                      "scipy",
                      "atlas_splitter>=0.1.4"
                      ],
    packages=find_packages(),
    scripts=[
        "bin/check-density-consistency",
        "bin/densities-across-regions",
        "bin/adjust-densities-across-regions",
        "bin/extract-mean-cell-counts"
    ],
    python_requires=">=3.6",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)

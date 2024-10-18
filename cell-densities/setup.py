#!/usr/bin/env python
# SPDX-License-Identifier: Apache-2.0

import imp
import sys

from setuptools import setup, find_packages

if sys.version_info < (3, 8):
    sys.exit("Sorry, Python < 3.6 is not supported")

VERSION = imp.load_source("", "nrrdhlp/version.py").__version__

setup(
    name="cell-density-validations",
    author=["Michael W. Reimann"],
    author_email="mwr@reimann.science",
    version=VERSION,
    description="validate cell densities",
    long_description="validate prescribed and built cell densities",
    url="http://bluebrain.epfl.ch",
    license="Apache-2.0",
    install_requires=["numpy >= 1.20.0",
                      "pandas >= 1.2.4",
                      "voxcell >= 3.1.5",
                      "scipy >= 1.8.0",
                      "atlas_splitter>=0.1.4"
                      ],
    packages=find_packages(),
    scripts=[
        "bin/check-density-consistency",
        "bin/densities-across-regions",
        "bin/adjust-densities-across-regions",
        "bin/extract-mean-cell-counts",
        "bin/extract-depth-profiles",
        "bin/validate-depth-profiles"
    ],
    python_requires=">=3.8",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)

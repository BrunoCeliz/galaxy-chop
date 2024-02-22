# This file is part of
# the galaxy-chop project (https://github.com/vcristiani/galaxy-chop)
# Copyright (c) Cristiani, et al. 2021, 2022, 2023
# License: MIT
# Full Text: https://github.com/vcristiani/galaxy-chop/blob/master/LICENSE.txt


# from distutils.core import setup
# from distutils.extension import Extension
# from Pyrex.Distutils import build_ext

"""C Implementations."""

from setuptools import Extension, find_packages, setup

setup(
    name="pot",
    version="0.1.0",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Natural Language :: English",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: Implementation :: CPython",
    ],
    ext_modules=[
        Extension(
            # the qualified name of the extension module to build
            "pot.pot",
            # the files to compile into our module relative to ``setup.py``
            ["potential_wrapped.c"],
        ),
    ],
)

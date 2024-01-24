# This file is part of
# the galaxy-chop project (https://github.com/vcristiani/galaxy-chop)
# Copyright (c) Cristiani, et al. 2021, 2022, 2023
# License: MIT
# Full Text: https://github.com/vcristiani/galaxy-chop/blob/master/LICENSE.txt

# =============================================================================
# DOCS
# =============================================================================

"""
GalaxyChop.

Implementation of a few galaxy dynamic decomposition methods.
"""

# =============================================================================
# IMPORTS
# =============================================================================

import importlib_metadata

NAME = "galaxychop"

VERSION = importlib_metadata.version(NAME)

__version__ = tuple(VERSION.split("."))

from . import constants, models, preproc, utils
from .core import (
    Galaxy,
    NoGravitationalPotentialError,
    ParticleSet,
    ParticleSetType,
    mkgalaxy,
)
from .io import read_hdf5, to_hdf5
from .pipeline import GchopPipeline, mkpipe


__all__ = [
    "Galaxy",
    "ParticleSet",
    "ParticleSetType",
    "NoGravitationalPotentialError",
    "io",
    "models",
    "preproc",
    "utils",
    "mkgalaxy",
    "read_hdf5",
    "to_hdf5",
    "constants",
    "GchopPipeline",
    "mkpipe",
]

del importlib_metadata

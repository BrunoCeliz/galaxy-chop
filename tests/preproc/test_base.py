# This file is part of
# the galaxy-chop project (https://github.com/vcristiani/galaxy-chop)
# Copyright (c) Cristiani, et al. 2021, 2022, 2023
# License: MIT
# Full Text: https://github.com/vcristiani/galaxy-chop/blob/master/LICENSE.txt

# WIP

import galaxychop as gchop

import numpy as np

import pandas as pd

import pytest


# =============================================================================
# TRANSFORMER ABC
# =============================================================================

# Bruno:
# Ojo cómo resulta el match entre lo del test_base de models/...
@pytest.mark.model
def GalaxyTransformerABC_not_implemethed():
    class Transformer(gchop.preproc.GalaxyTransformerABC):
        def transform(self, galaxy):
            return super().transform(galaxy)

        def checker(self, galaxy):
            return super().checker(galaxy)

    transformer = Transformer()

    with pytest.raises(NotImplementedError):
        transformer.transform(None)

    with pytest.raises(NotImplementedError):
        transformer.checker(None)


@pytest.mark.model
def test_GalaxyTransformerABC_repr():
    class Transformer(gchop.preproc.GalaxyTransformerABC):
        other = gchop.preproc.hparam(default=1)

        def transform(self, galaxy):
            ...

        def checker(self, galaxy):
            ...

    transformer = Transformer(other="zaraza")
    result = repr(transformer)
    expected = "Transformer(other='zaraza')"
    # Bruno: Ojo con esto...

    assert result == expected



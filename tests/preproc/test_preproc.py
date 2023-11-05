# This file is part of
# the galaxy-chop project (https://github.com/vcristiani/galaxy-chop)
# Copyright (c) Cristiani, et al. 2021, 2022, 2023
# License: MIT
# Full Text: https://github.com/vcristiani/galaxy-chop/blob/master/LICENSE.txt


# =============================================================================
# DOCS
# =============================================================================


"""Test utilities  galaxychop.preproc"""


# =============================================================================
# IMPORTS
# =============================================================================


from galaxychop import preproc

import pandas as pd


# =============================================================================
# TESTS
# =============================================================================


# Bruno:
# Aunque esto sea por completitud. Hay que deprecar estas funcs.
# *Como las funcs viejas siguen existiendo, no hace falta
# usar los transform ¿O sería lo mejor?
def test_center_and_align(galaxy):
    gal = galaxy(
        seed=42,
        stars_potential=True,
        dm_potential=True,
        gas_potential=True,
    )

    result = preproc.center_and_align(gal).to_dataframe()
    expected = preproc.salign.star_align(
        preproc.pcenter.center(gal)
    ).to_dataframe()

    pd.testing.assert_frame_equal(result, expected)


def test_is_centered_and_aligned(galaxy):
    gal = galaxy(
        seed=42,
        stars_potential=True,
        dm_potential=True,
        gas_potential=True,
    )

    assert preproc.is_centered_and_aligned(gal) is False
    assert (
        preproc.is_centered_and_aligned(preproc.pcenter.center(gal)) is False
    )
    assert (
        preproc.is_centered_and_aligned(preproc.salign.star_align(gal))
        is False
    )
    assert preproc.is_centered_and_aligned(
        preproc.salign.star_align(preproc.pcenter.center(gal))
    )
    assert preproc.is_centered_and_aligned(preproc.center_and_align(gal))

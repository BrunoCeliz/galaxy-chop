# This file is part of
# the galaxy-chop project (https://github.com/vcristiani/galaxy-chop)
# Copyright (c) Cristiani, et al. 2021, 2022, 2023
# License: MIT
# Full Text: https://github.com/vcristiani/galaxy-chop/blob/master/LICENSE.txt

# =============================================================================
# DOCS
# =============================================================================

"""Test utilities  galaxychop.preproc.smr_crop"""

# =============================================================================
# IMPORTS
# =============================================================================

from galaxychop.preproc import smr_crop

import numpy as np
import pytest


# =============================================================================
# TESTS
# =============================================================================


def test_half_star_mass_radius_crop(read_hdf5_galaxy):
    gal = read_hdf5_galaxy("gal394242.h5")

    cgal = smr_crop.half_star_mass_radius_crop(gal, num_radii=1)

    gal_m_s = gal.stars.m.sum().to_value()
    cgal_m_s = cgal.stars.m.sum().to_value()

    assert np.isclose(
        cgal_m_s, (gal_m_s / 2.0), rtol=1e-05, atol=1e-08, equal_nan=False
    )


# Bruno:
# Me "inspiro" de lo que hice en smr_crop.py
def test_gal_crop_invalid_num_radii(galaxy):
    gal = galaxy(seed=42)

    with pytest.raises(ValueError):
        smr_crop.half_star_mass_radius_crop(gal, num_radii=-1)


def test_is_gal_cropped_real_galaxy(read_hdf5_galaxy):
    gal = read_hdf5_galaxy("gal394242.h5")

    cgal = smr_crop.half_star_mass_radius_crop(gal)

    assert not smr_crop.is_star_cutted(gal, num_radii=1)
    assert smr_crop.is_star_cutted(cgal, num_radii=1)


def test_cutter_transformer(galaxy):
    gal = galaxy(seed=42)
    cutter = smr_crop.Cutter()

    assert cutter.transform(gal) == smr_crop.half_star_mass_radius_crop(gal)


def test_cutter_default_num_radii(galaxy):
    gal = galaxy("gal394242.h5")
    cutter = smr_crop.Cutter()

    assert cutter.transform(gal) == smr_crop.half_star_mass_radius_crop(
        gal, num_radii=3
    )


def test_cutter_notdefault_num_radii(galaxy):
    gal = galaxy("gal394242.h5")
    cutter = smr_crop.Cutter(num_radii=1)

    assert cutter.transform(gal) == smr_crop.half_star_mass_radius_crop(
        gal, num_radii=1
    )


def test_cutter_checker(galaxy):
    gal = galaxy(seed=42)
    cutter = smr_crop.Cutter()
    cgal = cutter.transform(gal)

    assert (cutter.checker(gal)) == (smr_crop.is_star_cutted(gal))
    assert (cutter.checker(cgal)) == (smr_crop.is_star_cutted(cgal))

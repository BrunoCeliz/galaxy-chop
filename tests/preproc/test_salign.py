# This file is part of
# the galaxy-chop project (https://github.com/vcristiani/galaxy-chop)
# Copyright (c) Cristiani, et al. 2021, 2022, 2023
# License: MIT
# Full Text: https://github.com/vcristiani/galaxy-chop/blob/master/LICENSE.txt

# =============================================================================
# DOCS
# =============================================================================

"""Test utilities  galaxychop.preproc.salign"""

# =============================================================================
# IMPORTS
# =============================================================================

from galaxychop.preproc import salign


import pytest


# =============================================================================
# ALIGN
# =============================================================================


def test_star_align_rcur0dot9(galaxy):
    gal = galaxy(seed=42)

    agal = salign.star_align(gal, r_cut=0.9)

    df = gal.to_dataframe()
    adf = agal.to_dataframe()

    changed = [
        "x",
        "y",
        "z",
        "vx",
        "vy",
        "vz",
        "Jx",
        "Jy",
        "Jz",
        "kinetic_energy",
        "total_energy",
    ]

    for colname in df.columns[~df.columns.isin(changed)]:
        ocol = df[colname]
        acol = adf[colname]
        assert (ocol == acol).all(), colname

    for colname in changed:
        ocol = df[colname]
        acol = adf[colname]
        assert not (ocol == acol).all(), colname


def test_star_align(galaxy):
    gal = galaxy(seed=42)

    agal = salign.star_align(gal)

    df = gal.to_dataframe()
    adf = agal.to_dataframe()

    changed = [
        "x",
        "y",
        "z",
        "vx",
        "vy",
        "vz",
        "Jx",
        "Jy",
        "Jz",
        "kinetic_energy",
        "total_energy",
    ]

    for colname in df.columns[~df.columns.isin(changed)]:
        ocol = df[colname]
        acol = adf[colname]
        assert (ocol == acol).all(), colname

    for colname in changed:
        ocol = df[colname]
        acol = adf[colname]
        assert not (ocol == acol).all(), colname


def test_star_align_invalid_rcut(galaxy):
    gal = galaxy(seed=42)

    with pytest.raises(ValueError):
        salign.star_align(gal, r_cut=-1)


def test_is_star_aligned_real_galaxy(read_hdf5_galaxy):
    gal = read_hdf5_galaxy("gal394242.h5")

    agal = salign.star_align(gal, r_cut=5)

    assert not salign.is_star_aligned(gal, r_cut=5)
    assert salign.is_star_aligned(agal, r_cut=5)


def test_is_star_aligned_fake_galaxy(galaxy):
    gal = galaxy(seed=42)

    agal = salign.star_align(gal, r_cut=5)

    assert not salign.is_star_aligned(gal, r_cut=5)
    assert salign.is_star_aligned(agal, r_cut=5)


# Bruno:
# Las funcs privs no se testean...
# Para testear las clases, las inicializamos y comparamos
# que sus métodos hagan lo mismo que las funciones aparte
def test_aligner_transformer(galaxy):
    gal = galaxy(seed=42)
    aligner = salign.Aligner()

    # Bruno: Ídem que pcenter...
    class_agal = aligner.transform(gal)
    class_df = class_agal.to_dataframe()

    func_agal = salign.star_align(gal)
    func_df = func_agal.to_dataframe()

    assert class_df.equals(func_df)


# Bruno:
# Añado un test para saber si inicializa bien el r_cut:
def test_aligner_default_r_cut(read_hdf5_galaxy):
    gal = read_hdf5_galaxy("gal394242.h5")
    aligner = salign.Aligner()

    class_agal = aligner.transform(gal)
    class_df = class_agal.to_dataframe()

    func_agal = salign.star_align(gal, r_cut=30)
    func_df = func_agal.to_dataframe()

    assert class_df.equals(func_df)


# Bruno:
# Y otro para saber si pone bien un r_cut nuevo:
def test_aligner_notdefault_r_cut(read_hdf5_galaxy):
    gal = read_hdf5_galaxy("gal394242.h5")
    aligner = salign.Aligner(r_cut=10)

    class_agal = aligner.transform(gal)
    class_df = class_agal.to_dataframe()

    func_agal = salign.star_align(gal, r_cut=10)
    func_df = func_agal.to_dataframe()

    assert class_df.equals(func_df)


def test_aligner_checker(galaxy):
    gal = galaxy(seed=42)
    aligner = salign.Aligner()
    agal = aligner.transform(gal)

    assert (aligner.checker(gal)) == (salign.is_star_aligned(gal))
    assert (aligner.checker(agal)) == (salign.is_star_aligned(agal))

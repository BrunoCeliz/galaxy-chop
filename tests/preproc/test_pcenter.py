# This file is part of
# the galaxy-chop project (https://github.com/vcristiani/galaxy-chop)
# Copyright (c) Cristiani, et al. 2021, 2022, 2023
# License: MIT
# Full Text: https://github.com/vcristiani/galaxy-chop/blob/master/LICENSE.txt


# =============================================================================
# DOCS
# =============================================================================


"""Test utilities  galaxychop.preproc.pcenter"""


# =============================================================================
# IMPORTS
# =============================================================================


from galaxychop.preproc import pcenter

import pytest


# =============================================================================
# CENTER
# =============================================================================


@pytest.mark.parametrize("with_potential", [True, False])
def test_center_without_potential_energy(galaxy, with_potential):
    gal = galaxy(
        seed=42,
        stars_potential=False,
        dm_potential=False,
        gas_potential=False,
    )
    if with_potential:
        with pytest.raises(ValueError):
            pcenter.center(gal)
    else:
        cgal = pcenter.center(gal)

        df = gal.to_dataframe()
        cdf = cgal.to_dataframe()

        changed = ["x", "y", "z", "vx", "vy", "vz", "Jx", "Jy", "Jz"]

        for colname in df.columns[~df.columns.isin(changed)]:
            ocol = df[colname]
            ccol = cdf[colname]
            assert (ocol == ccol).all()

        for colname in changed:
            ocol = df[colname]
            ccol = cdf[colname]
            assert not (ocol == ccol).all()


def test_center(galaxy):
    gal = galaxy(
        seed=42,
        stars_potential=True,
        dm_potential=True,
        gas_potential=True,
    )

    cgal = pcenter.center(gal)

    df = gal.to_dataframe()
    cdf = cgal.to_dataframe()

    changed = ["x", "y", "z", "vx", "vy", "vz", "Jx", "Jy", "Jz"]

    for colname in df.columns[~df.columns.isin(changed)]:
        ocol = df[colname]
        ccol = cdf[colname]
        assert (ocol == ccol).all()

    for colname in changed:
        ocol = df[colname]
        ccol = cdf[colname]
        assert not (ocol == ccol).all()


def test_is_centered_without_potential_energy(galaxy):
    gal = galaxy(
        seed=42,
        stars_potential=False,
        dm_potential=False,
        gas_potential=False,
    )
    with pytest.raises(ValueError):
        pcenter.is_centered(gal)


def test_is_centered(galaxy):
    gal = galaxy(
        seed=42,
        stars_potential=True,
        dm_potential=True,
        gas_potential=True,
    )

    cgal = pcenter.center(gal)

    assert not pcenter.is_centered(gal)
    assert pcenter.is_centered(cgal)


# Bruno:
# Para testear las clases, las inicializamos y comparamos
# que sus métodos hagan lo mismo que las funciones aparte
# ¿Como argumento tengo que poner a la clase?
def test_centralizer_transformer(galaxy):
    gal = galaxy(
        seed=42,
        stars_potential=True,
        dm_potential=True,
        gas_potential=True,
    )

    centralizer = pcenter.Centralizer()
    class_cgal = centralizer.transform(gal)
    class_df = class_cgal.to_dataframe()

    func_cgal = pcenter.center(gal)
    func_df = func_cgal.to_dataframe()

    assert class_df.equals(func_df)


def test_centralizer_checker(galaxy):
    gal = galaxy(
        seed=42,
        stars_potential=True,
        dm_potential=True,
        gas_potential=True,
    )

    centralizer = pcenter.Centralizer()
    cgal = centralizer.transform(gal)

    assert (centralizer.checker(gal)) == (pcenter.is_centered(gal))
    assert (centralizer.checker(cgal)) == (pcenter.is_centered(cgal))

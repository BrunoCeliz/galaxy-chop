# This file is part of
# the galaxy-chop project (https://github.com/vcristiani/galaxy-chop)
# Copyright (c) Cristiani, et al. 2021, 2022, 2023
# License: MIT
# Full Text: https://github.com/vcristiani/galaxy-chop/blob/master/LICENSE.txt

# WIP (¿Hacer folders para C y fortran?)

# =============================================================================
# DOCS
# =============================================================================

"""Test utilities  galaxychop.preproc.grispy_potential"""

# =============================================================================
# IMPORTS
# =============================================================================

from galaxychop.preproc import grispy_potential


import numpy as np

# =============================================================================
# TESTS
# =============================================================================

# Bruno:
# Probemos: i) si l_box (tamaño del grid) encierra a todas las partículas; \
# ii) si el grid generado es un objeto de GriSPy; iii) (siguiendo los métodos \
# de otros test, revisar) si el return de la func es un ndarray de \
# shape = (n,1); ¿Algo más? ¿O como es un import de otra librería no hace
# falta hacer tests grained?


def test_make_grid(galaxy):
    gal = galaxy(
        seed=42,
        stars_potential=False,
        dm_potential=False,
        gas_potential=False,
    )

    df = gal.to_dataframe()

    # ¿Así? ¿O mejor me creo arrays de shape = (n,1)?
    x_sys = df.x.values
    y_sys = df.y.values
    z_sys = df.z.values

    l_box, grid = grispy_potential.make_grid(x_sys, y_sys, z_sys)

    # Bruno:
    # Probamos...
    assert np.max(x_sys) < l_box
    assert np.min(x_sys) > -l_box
    assert np.max(y_sys) < l_box
    assert np.min(y_sys) > -l_box
    assert np.max(z_sys) < l_box
    assert np.min(z_sys) > -l_box

    # Para lo del grid, que la cantidad de pts en el grid sea igual \
    # a la cantidad de partículas en el sistema (creo que es un \
    # buen sanity-check);
    # *Existe el "grid.contains(particle)", pero no lo uso acá
    assert len(x_sys) == grid.ndata


def test_potential_grispy(galaxy):
    gal = galaxy(
        seed=42,
        stars_potential=False,
        dm_potential=False,
        gas_potential=False,
    )

    df = gal.to_dataframe()

    x_sys = df.x.values
    y_sys = df.y.values
    z_sys = df.z.values
    m_sys = df.m.values

    l_box, grid = grispy_potential.make_grid(x_sys, y_sys, z_sys)

    # Bruno:
    # Ojo que hasta acá repito todo lo de la func anterior...

    centre = np.array([x_sys[0], y_sys[0], z_sys[0]])

    softening = 0.5

    epot = grispy_potential.potential_grispy(
        centre,
        x_sys,
        y_sys,
        z_sys,
        m_sys,
        softening,
        5 * softening,
        0.1 * l_box,
        l_box,
        grid,
    )

    # Bruno:
    # Probamos...
    assert len(epot) == 1
    assert epot < 0

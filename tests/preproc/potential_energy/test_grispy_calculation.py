# This file is part of
# the galaxy-chop project (https://github.com/vcristiani/galaxy-chop)
# Copyright (c) Cristiani, et al. 2021, 2022, 2023
# License: MIT
# Full Text: https://github.com/vcristiani/galaxy-chop/blob/master/LICENSE.txt


# =============================================================================
# DOCS
# =============================================================================


"""Test utilities  galaxychop.preproc.potential_energy.grispy_calculation"""


# =============================================================================
# IMPORTS
# =============================================================================


from galaxychop.preproc.potential_energy import grispy_calculation

import grispy as gsp

import numpy as np


# =============================================================================
# TESTS
# =============================================================================

# Bruno:
# Probemos:
# i) si l_box (tamaño del grid) encierra a todas las partículas;
# ii) si el grid generado es un objeto de GriSPy;
# iii) (revisar los métodos de otros test) si el return de la func es un
# ndarray de shape = (n,1);
# ¿Algo más? ¿O como es un import de otra librería no hace falta hacer tests
# más "grained"?


# Bruno:
# Pruebo de otra forma para usar glxs reales...
# def test_make_grid(galaxy):
#     gal = galaxy(
#         seed=42,
#         stars_potential=False,
#         dm_potential=False,
#         gas_potential=False,
#     )
def test_make_grid(read_hdf5_galaxy):
    gal = read_hdf5_galaxy("gal394242.h5")
    df = gal.to_dataframe()

    # ¿Así? ¿O mejor me creo arrays de shape = (n,1)?
    x_sys = df.x.values
    y_sys = df.y.values
    z_sys = df.z.values

    l_box, grid = grispy_calculation.make_grid(x_sys, y_sys, z_sys)

    # Bruno: Que estén todas las partículas encerradas
    # ¿Debería usar más los "and"?
    assert l_box > 0
    assert np.max(x_sys) < l_box
    assert np.min(x_sys) > -l_box
    assert np.max(y_sys) < l_box
    assert np.min(y_sys) > -l_box
    assert np.max(z_sys) < l_box
    assert np.min(z_sys) > -l_box

    # Para lo del grid, que la cantidad de pts en el grid sea igual
    # a la cantidad de partículas en el sistema (creo que es un
    # buen sanity-check);
    # *Existe el "grid.contains(particle)", pero no lo uso acá
    assert len(x_sys) == grid.ndata
    # Bruno: Ojo...
    assert isinstance(grid, gsp.GriSPy)


def test_potential_grispy_one_particle(read_hdf5_galaxy):
    gal = read_hdf5_galaxy("gal394242.h5")

    df = gal.to_dataframe()

    x_sys = df.x.values
    y_sys = df.y.values
    z_sys = df.z.values
    m_sys = df.m.values

    l_box, grid = grispy_calculation.make_grid(x_sys, y_sys, z_sys)

    # Bruno: Mismo procedimiento para armar el grid, pero testeo
    # que funcione para 1 partícula (teniendo en cuenta el resto)
    centre = np.column_stack((x_sys[0], y_sys[0], z_sys[0]))
    softening = 0.5

    epot = grispy_calculation.potential_grispy(
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

    # Bruno: ¿Suficiente? Porque la "calidad" del cómputo
    # lo probamos en el init de potential_energy...
    assert isinstance(epot, float)
    assert epot < 0

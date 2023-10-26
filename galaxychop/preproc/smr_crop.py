# This file is part of
# the galaxy-chop project (https://github.com/vcristiani/galaxy-chop)
# Copyright (c) Cristiani, et al. 2021, 2022, 2023
# License: MIT
# Full Text: https://github.com/vcristiani/galaxy-chop/blob/master/LICENSE.txt

# =============================================================================
# DOCS
# =============================================================================

"""Utilities for crop a galaxy based on stellar mass."""

# =============================================================================
# IMPORTS
# =============================================================================

import numpy as np

from ..core import data

# =============================================================================
# API
# =============================================================================


# Bruno:
# "smr" == "Stellar Mass Radius". Esta función no debería ser privada, ya que
# esta forma de medir el tamaño de la galaxia es re útil y el usuario puede
# requerirla en caso de que no tenga el dato de la simulación.
# Dicho sea de paso, exsiten otras formas de medir la galaxia e.g. "Half Mass
# Radius" (considerando TODAS las partículas (gas, DM y stars)); y más aún,
# también le podemos devolver el valor de la MASA (no sólo el radio)
# encerrada, como para evitar que el usuario escriba esas líneas de código
# aparte, sin usar GlxChop.
# Otra cosa, aunque quizá sea innecesario, debería existir una columna de
# "r" para las partículas: si bien ya existe "x, y, z", usamos para muchas
# cosas la distancia galactocéntrica "r" que debería incorporarse a la clase
# Galaxy (o dentro del mismo ParticleSet) (!).
def _get_half_smr_crop(sdf, cut_radius_factor):
    sdf = sdf[["x", "y", "z", "m"]].copy()

    sdf["radius"] = np.sqrt(sdf.x**2 + sdf.y**2 + sdf.z**2)
    sdf.drop(["x", "y", "z"], axis="columns", inplace=True)

    sdf.sort_values("radius", inplace=True)
    # Bruno:
    # Por esto mismo lo digo, acá ya calculaste el "r", ¿Por qué no
    # agregárselo a la Galaxia?

    sdf["m_cumsum"] = sdf.m.cumsum()
    # Más aún, acá calculaste la masa -> podría guardarse y devolverla
    # como dato de la galaxia al usuario...
    sdf.drop(["m"], axis="columns", inplace=True)

    half_m_cumsum = sdf.iloc[-1].m_cumsum / 2
    sdf["half_m_cumsum_diff"] = np.abs(sdf.m_cumsum - half_m_cumsum)

    cut_radius = (
        sdf.iloc[sdf.half_m_cumsum_diff.argmin()].radius * cut_radius_factor
    )

    cut_idxs = sdf[sdf.radius > cut_radius].index.to_numpy()

    del sdf

    return cut_idxs, cut_radius


def half_star_mass_radius_crop(galaxy, *, num_radii=3):
    """
    Crop select stars within a specified number of the radii enclosing \
    half fractions of the stellar mass.

    Parameters
    ----------
    galaxy : galaxychop.Galaxy
        The galaxy object for which to determine half of the
        mass-enclosing radii.
    num_radii : int, optional
        The number of radii to consider. Default is 3.

    Returns
    -------
    galaxychop.Galaxy
        A new galaxy object containing stars within the specified radii
        enclosing various fractions of the stellar mass.

    """
    # We convert the stars into a dataframe
    stars_df = galaxy.stars.to_dataframe()

    # We check which rows to delete and what cutoff radius it gives us
    to_trim_idxs, _ = _get_half_smr_crop(stars_df, cut_radius_factor=num_radii)
    trim_stars_df = stars_df.drop(to_trim_idxs, axis="rows")

    # We create a new particle set with the new stars.
    trim_stars = data.ParticleSet(
        ptype=data.ParticleSetType.STARS,
        m=trim_stars_df["m"].to_numpy(),
        x=trim_stars_df["x"].to_numpy(),
        y=trim_stars_df["y"].to_numpy(),
        z=trim_stars_df["z"].to_numpy(),
        vx=trim_stars_df["vx"].to_numpy(),
        vy=trim_stars_df["vy"].to_numpy(),
        vz=trim_stars_df["vz"].to_numpy(),
        potential=trim_stars_df["potential"].to_numpy(),
        softening=galaxy.stars.softening,
    )

    del trim_stars_df

    # Bruno:
    # Si cortamos las estrellas hasta 3r_half, ¿Por qué no al gas? ;
    # Además, también se podrían agregar otros cortes (e.g. Einasto, \
    # NFW, half light radii (!, aunque para ello debería considerarse \
    # la magnitud/luminosidad de cada partícula, que si bien no es \
    # imposible pero estira el scope de esto y habría que retocar \
    # muchas otras varias cosas...))
    dm = galaxy.dark_matter.copy()
    gas = galaxy.gas.copy()

    # D: aca es cuando hace una copia/ crea una galaxia no?

    trim_galaxy = data.Galaxy(stars=trim_stars, dark_matter=dm, gas=gas)

    return trim_galaxy

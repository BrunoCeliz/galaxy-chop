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

from ._base import GalaxyTransformerABC, hparam
from ..core import data
from ..utils import doc_inherit

# =============================================================================
# INTERNALS
# =============================================================================
# Bruno: ¿?; btw: esta func debería ser pública/calcular métricas de la galaxia
# y anexarla a la misma...


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


# =============================================================================
# CUTTER CLASS
# =============================================================================
# Bruno:
# "Cropper" ~> granjero; "Chopper" ~> Helicóptero o Moto;
# "Cutter" ~> re de nene; lpm


class Cutter(GalaxyTransformerABC):
    """
    Cutter class.

    Given the positions and mass of particles, compute the stellar
    half mass radius (radii of the sphere centered at the origin that
    encloses half of the total sum of stellar particles mass) and
    return only the stellar component inside of a multiple of it.

    """

    num_radii = hparam(default=3)
    # Bruno: Nos suelen gustar ~30 kpc. Pero si 3 r_half es mucho menor,
    # hay que tener cuidado...

    @doc_inherit(GalaxyTransformerABC.transform)
    def transform(self, galaxy, num_radii):
        # Bruno:
        # Esto returnea glx, pero el dato del r_half es más importante
        return half_star_mass_radius_crop(galaxy, num_radii)

    # Bruno:
    # No un "checker" como tal, habría que hacerlo. Un buen test es que
    # la distancia galactocéntica máxima no sea mayor a
    # num_radii*r_half...
    @doc_inherit(GalaxyTransformerABC.checker)
    def checker(self, galaxy, **kwargs):
        # Bruno:
        # ¿Mejor forma que no sea con el "**"? ¿Aclarar? Ojo...
        return is_star_cutted(galaxy, **kwargs)


# =============================================================================
# API FUNCTIONS
# =============================================================================
# Bruno: Rev como dejar este título de sec como en otros módulos...


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


def is_star_cutted(galaxy, *, num_radii=3, rtol=1e-05, atol=1e-08):
    """
    Validate if the galaxy is cutted.

    Parameters
    ----------
    galaxy : ``Galaxy class`` object
    num_radii : float, optional
        Default value =  3. If it's provided,
        it must be positive and the galaxy is
        cutted at that radii, with only particles
        at radii smaller than num_radii*r_half.
    rtol : float
        Relative tolerance. Default value = 1e-05.
    atol : float
        Absolute tolerance. Default value = 1e-08.

    Returns
    -------
    bool
        True if galaxy is cutted at given radii, and
        there is no stellar particles at further distances,
        False otherwise.

    """
    # Now we extract only the needed column to rotate the galaxy
    stars_df = galaxy.stars.to_dataframe(attributes=["m", "x", "y", "z"])

    to_trim_idxs, cut_radius = _get_half_smr_crop(stars_df, num_radii)
    trim_stars_df = stars_df.drop(to_trim_idxs, axis="rows")

    # Distances of all stellar particles that "survives"
    distances = np.sqrt(
        trim_stars_df.x**2 + trim_stars_df.y**2 + trim_stars_df.z**2
    )

    # maximum distance index of all particles
    maxdist_idx = np.argmin(distances)[0]
    max_values = distances[maxdist_idx]

    # Bruno:
    # Probando cómo checkear ~> todas las no cortadas no deben estar
    # a mayor distancia que el corte (y como test también agregar que
    # las cortadas deben estar más lejos que el corte (!!!))
    return np.all(np.less([max_values, cut_radius]))

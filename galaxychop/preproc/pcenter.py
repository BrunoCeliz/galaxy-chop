# This file is part of
# the galaxy-chop project (https://github.com/vcristiani/galaxy-chop)
# Copyright (c) Cristiani, et al. 2021, 2022, 2023
# License: MIT
# Full Text: https://github.com/vcristiani/galaxy-chop/blob/master/LICENSE.txt

# =============================================================================
# DOCS
# =============================================================================

"""Utilities to center a galaxy."""

# Bruno:
# Falta terminar la documentación del nuevo center (!!!)

# =============================================================================
# IMPORTS
# =============================================================================

import numpy as np

from ._base import GalaxyTransformerABC
from ..core import data
from ..utils import doc_inherit

# =============================================================================
# CENTRALIZER CLASS
# =============================================================================


# Bruno:
# Agregar que acomode las velocidades (no es obvio, pero si está centrado
# a datos de simulaciones es obligatorio que el v_cm = 0 por los
# cálculos de Jx, Jym Jz...) (!!!)
class Centralizer(GalaxyTransformerABC):
    # Bruno:
    # Ojo con la doc...
    """
    Centralizer class.

    Given the positions and potential energy of particles, check and
    center their positions relative to the system.

    """

    @doc_inherit(GalaxyTransformerABC.transform)
    def transform(self, galaxy, **kwargs):
        return center(galaxy, **kwargs)

    @doc_inherit(GalaxyTransformerABC.checker)
    def checker(self, galaxy, **kwargs):
        return is_centered(galaxy, **kwargs)


# =============================================================================
# API FUNCTIONS
# =============================================================================


def center(galaxy, with_potential=True):
    """
    Galaxy particle centering.

    Centers the position of all galaxy particles respect to the position of the
    lowest potential particle. If with_potential is True, it subtracts the
    position columns of the galaxy dataframe by the position with the lowest
    potential value. If with_potential is False, it subtracts the position
    columns of the galaxy dataframe by the geometric center of the positions.

    Parameters
    ----------
    galaxy : ``Galaxy class`` object

    with_potential : boolean. True by default

    Returns
    -------
    galaxy : new ``Galaxy class`` object
        A new galaxy object with centered positions respect to the position of
        the lowest potential particle.

    """

    # Bruno:
    # Agrego la opción de centrar según centroide =/= partícula de min pot
    # TO DO: Que acomode las velocidades (al mismo tiempo o por separado,
    # ver cómo implementar...)
    if not galaxy.has_potential_ and with_potential:
        raise ValueError(
            "Galaxy must has the potential energy. Otherwise, use \
            with_potential = False"
        )

    if with_potential:
        # We extract only the needed column to centrer the galaxy
        df = galaxy.to_dataframe(
            attributes=["ptypev", "x", "y", "z", "vx", "vy", "vz", "potential"]
        )

        # minimum potential index of all particles and we extract data
        # frame row
        minpot_idx = df.potential.argmin()
        min_values = df.iloc[minpot_idx]

        # We subtract all position columns by the position with the lowest
        # potential value and replace this new position columns on dataframe
        columns = ["x", "y", "z"]
        df[columns] = df[columns] - min_values[columns]

    else:
        # We use only positions
        df = galaxy.to_dataframe(
            attributes=["ptypev", "x", "y", "z", "vx", "vy", "vz"]
        )

        # Compute the geometric center
        x_cm = np.mean(df["x"].values)
        y_cm = np.mean(df["y"].values)
        z_cm = np.mean(df["z"].values)

        # We subtract all position columns by the new origin
        # and replace on dataframe
        df.loc[:, "x"] -= x_cm
        df.loc[:, "y"] -= y_cm
        df.loc[:, "z"] -= z_cm

    # Bruno:
    # Nuevo -> Acomodo según el v_CM (con promedios en realidad...)

    # Compute the velocity of the galaxy within the cosmological box
    vx_cm = np.mean(df["vx"].values)
    vy_cm = np.mean(df["vy"].values)
    vz_cm = np.mean(df["vz"].values)

    # And modify the dataframe
    df.loc[:, "vx"] -= vx_cm
    df.loc[:, "vy"] -= vy_cm
    df.loc[:, "vz"] -= vz_cm

    # We split the dataframe by particle type.
    stars = df[df.ptypev == data.ParticleSetType.STARS.value]
    dark_matter = df[df.ptypev == data.ParticleSetType.DARK_MATTER.value]
    gas = df[df.ptypev == data.ParticleSetType.GAS.value]

    # patch
    new = galaxy.disassemble()

    new.update(
        x_s=stars.x.to_numpy(),
        y_s=stars.y.to_numpy(),
        z_s=stars.z.to_numpy(),
        x_dm=dark_matter.x.to_numpy(),
        y_dm=dark_matter.y.to_numpy(),
        z_dm=dark_matter.z.to_numpy(),
        x_g=gas.x.to_numpy(),
        y_g=gas.y.to_numpy(),
        z_g=gas.z.to_numpy(),
    )

    return data.mkgalaxy(**new)


# Bruno:
# Claro, ahora el chekcer también debería diferenciar si
# tiene o no tiene el potencial... Vamos a patearlo por ahora
# porque en una de esas no vale la pena considerar el caso de
# Galaxia sin potencial...
def is_centered(galaxy, *, rtol=1e-05, atol=1e-08):
    """
    Validate if the galaxy is centered.

    Parameters
    ----------
    galaxy : ``Galaxy class`` object
    rtol : float
        Relative tolerance. Default value = 1e-05.
    atol : float
        Absolute tolerance. Default value = 1e-08.

    Returns
    -------
    bool
        True if galaxy is centered respect to the position of the lowest
        potential particle, False otherwise.

    """
    if not galaxy.has_potential_:
        raise ValueError("Galaxy must has the potential energy.")

    # We extract only the needed column to center the galaxy
    # Bruno: ¡Agrego lo de las velocidades!
    df = galaxy.to_dataframe(
        attributes=["x", "y", "z", "vx", "vy", "vz", "potential"]
    )

    # minimum potential index of all particles and we extract data frame row
    minpot_idx = df.potential.argmin()
    min_values = df.iloc[minpot_idx]

    # mean velocity of the particles
    vx_cm = np.mean(df["vx"].values)
    vy_cm = np.mean(df["vy"].values)
    vz_cm = np.mean(df["vz"].values)

    # Bruno: Ojo... Testear que funque
    return (
        np.allclose(min_values[["x", "y", "z"]], 0, rtol=rtol, atol=atol)
    ) and (
        np.isclose(
            np.sqrt(vx_cm**2 + vy_cm**2 + vz_cm**2),
            0,
            rtol=rtol,
            atol=atol,
        )
    )

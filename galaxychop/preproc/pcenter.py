# This file is part of
# the galaxy-chop project (https://github.com/vcristiani/galaxy-chop)
# Copyright (c) Cristiani, et al. 2021, 2022, 2023
# License: MIT
# Full Text: https://github.com/vcristiani/galaxy-chop/blob/master/LICENSE.txt

# =============================================================================
# DOCS
# =============================================================================

"""Utilities to center a galaxy."""

# =============================================================================
# IMPORTS
# =============================================================================

import numpy as np

# Bruno:
# Me falta importar de las otras carpetas...
from ._base import GalaxyTransformerABC  # , hparam -> # Bruno: Unused (yet)
from ..core import data
from ..utils import doc_inherit

# =============================================================================
# CENTRALIZER CLASS
# =============================================================================


# Bruno:
# Algo importantísimo que me parece que falta es que, si la galaxia NO \
# tiene calculado el potencial, que centre según la media de las posiciones \
# en cada eje (y, calro, que se lo haga saber al usuario) i.e. si no quiero \
# calcular el potencial y soy confianzudo, hacé x_i = x_i - prom(x) para \
# todas las partículas...
# Más aún ¡Que también haga lo mismo para las velocidades! (no es obvio, \
# pero si está centrado a datos de simulaciones es obligatorio que el \
# v_cm = 0 por los cálculos de Jx, Jym Jz...) (!!!)
class Centralizer(GalaxyTransformerABC):
    # Bruno:
    # Ojo con la doc. NO sé cuál es el resultado final
    # (y/o qué correspondería a la clase/método (!)).
    """
    Centralizer class.

    Given the positions and potential energy of particles, check and
    center their positions relative to the system.

    """

    # Bruno:
    # Este método sólo llama a la func (priv) externa
    @doc_inherit(GalaxyTransformerABC.transform)
    def transform(self, galaxy):
        # Bruno: Cosa -> Así como para el cálculo de potencial,
        # esto "desarma" galaxias y manipula sus atributos...
        # En una pipeline no queremos ese behaviour (futurísimo work)...
        return center(galaxy)  # Así no más?

    # Bruno:
    # ¿Esta func la quiero como método? Como es un "centralizer",
    # lo dejemos adentro...
    @doc_inherit(GalaxyTransformerABC.checker)
    def checker(self, galaxy, **kwargs):
        # Bruno:
        # Again, quizás como warning, pero que no
        # corte la línea de producción...
        # (el cómputo del potencial suele molestar mucho para una muestra \
        # grande de galaxias...)
        # El asterisco acá adentro es lo que le molesta...
        # ¿Mejor forma que no sea con el "**"? ¿Aclarar? Ojo...
        return is_centered(galaxy, **kwargs)


# =============================================================================
# API FUNCTIONS
# =============================================================================
# Bruno: Rev como dejar este título de sec como en otros módulos...


def center(galaxy):
    """
    Galaxy particle centering.

    Centers the position of all galaxy particles respect to the position of the
    lowest potential particle.

    Parameters
    ----------
    galaxy : ``Galaxy class`` object

    Returns
    -------
    galaxy : new ``Galaxy class`` object
        A new galaxy object with centered positions respect to the position of
        the lowest potential particle.

    """
    if not galaxy.has_potential_:
        raise ValueError("galaxy must has the potential energy")

    # We extract only the needed column to centrer the galaxy
    df = galaxy.to_dataframe(attributes=["ptypev", "x", "y", "z", "potential"])

    # minimum potential index of all particles and we extract data frame row
    minpot_idx = df.potential.argmin()
    min_values = df.iloc[minpot_idx]

    # We subtract all position columns by the position with the lowest
    # potential value and replace this new position columns on dataframe
    columns = ["x", "y", "z"]
    df[columns] = df[columns] - min_values[columns]

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
    # ¿Por qué los "**" antes del "new"?
    return data.mkgalaxy(**new)


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
        raise ValueError("galaxy must has the potential energy")

    # We extract only the needed column to centrer the galaxy
    df = galaxy.to_dataframe(attributes=["x", "y", "z", "potential"])

    # minimum potential index of all particles and we extract data frame row
    minpot_idx = df.potential.argmin()
    min_values = df.iloc[minpot_idx]

    return np.allclose(min_values[["x", "y", "z"]], 0, rtol=rtol, atol=atol)

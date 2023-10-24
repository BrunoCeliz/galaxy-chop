# This file is part of
# the galaxy-chop project (https://github.com/vcristiani/galaxy-chop)
# Copyright (c) Cristiani, et al. 2021, 2022, 2023
# License: MIT
# Full Text: https://github.com/vcristiani/galaxy-chop/blob/master/LICENSE.txt

# =============================================================================
# DOCS
# =============================================================================

"""Preprocessing module."""

# =============================================================================
# IMPORTS
# =============================================================================

# Bruno:
# Cuando agregue (nuevamente) el Octree de C, volver a acá...
"""from ._base import (
    _PTYPES_ORDER,
    hparam,
    GalaxyTransformerABC,
)"""  # Bruno: ¿"imported but unused"? ¿Y el resto de imports?
from .pcenter import Centralizer
from .salign import Aligner
from .smr_crop import half_star_mass_radius_crop # D: Aun no lo hicimos clase
from .potential_energy import Potential # D: aca no se si esta bien llamado

"""from .potential_energy import (
    POTENTIAL_BACKENDS,
    DEFAULT_POTENTIAL_BACKEND,
    Potential,
    fortran_potential,
    grispy_potential,
    numpy_potential,
)"""  # Bruno: ¿"imported but unused"? ¿Y el resto de imports?

# Bruno:
# Acomodar...
# D: Esto como quedaria? que estructura se le da
__all__ = [
    # pcenter
    "center",
    "is_centered",
    "potential",
    "star_align",
    "is_star_aligned",
    "center_and_align",
    "half_star_mass_radius_crop",
]

# =============================================================================
# FUNCTIONS
# =============================================================================

# Bruno:
# ¿Esto era lo otro que molestaba? Deberíamos repetir un _base.py con un ABC \
# para estos pre-procesadores y unificar. Además, agregar el Octree + la \
# integración de C (Cython u otro)...
# ¡El r_cut puede ser un hparam!


def center_and_align(galaxy, *, r_cut=None):
    """
    Sequentially performs centering and alignment.

    ``center_and_align(galaxy) <==> star_align(center(galaxy))``

    Parameters
    ----------
    galaxy : ``Galaxy class`` object
    r_cut : float, optional
        Default value =  None. If it's provided, it must be positive and the
        rotation matrix `A` is calculated from the particles with radii smaller
        than r_cut.

    Returns
    -------
    galaxy: new ``Galaxy class`` object
        A new galaxy object with centered positions respect to the position of
        the lowest potential particle and their total angular momentum aligned
        with the z-axis.

    """
    # Bruno:
    # Remember que ahora son clases (!)
    center = Centralizer.transform()
    align = Aligner.transform()
    # aligned = Aling.transform(centered, r_cut=r_cut) # D: esto tendria que ir asi segun entiendo
    # Bruno: 1ro se inicializa el tranform, y luego se le da de comer.
    # Quizás haya otra forma, por ahora piso la galaxia que se come la func.
    galaxy = center(galaxy)  # Centering
    galaxy = align(galaxy, r_cut)  # Aligning

    return galaxy


# Bruno:
# Ídem para "rtol" y "atol"...
def is_centered_and_aligned(galaxy, *, r_cut=None, rtol=1e-05, atol=1e-08):
    """
    Validate if the galaxy is centered and aligned.

    ``is_center_and_align(galaxy) <==> is_centered(galaxy) and \
                                       is_star_aligned(galaxy)``

    Parameters
    ----------
    galaxy : ``Galaxy class`` object
    r_cut : float, optional
        Default value =  None. If it's provided, it must be positive and the
        rotation matrix `A` is calculated from the particles with radii smaller
        than r_cut.
    rtol : float
        Relative tolerance. Default value = 1e-05.
    atol : float
        Absolute tolerance. Default value = 1e-08.

    Returns
    -------
    bool
        True if galaxy is centered respect to the position of the lowest
        potential particle, and if the total angular momentum of the galaxy
        is aligned with the z-axis, False otherwise.

    """
    # Bruno:
    # Remember que ahora son clases (!)
    check_center = Centralizer.is_centered()
    check_align = Aligner.is_aligned()
    return check_center(
        galaxy, rtol=rtol, atol=atol) and check_align(
        galaxy, r_cut=r_cut, rtol=rtol, atol=atol
    )

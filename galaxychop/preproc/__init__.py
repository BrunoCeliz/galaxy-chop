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
# Acomodar todo con lo nuevo (Clases que llaman a func externas)...
from .pcenter import Centralizer, center, is_centered
from .potential_energy import (
    Potentializer,
    potential,
)  # D: aca no se si esta bien llamado
from .salign import Aligner, is_star_aligned, star_align
from .smr_crop import Cutter, half_star_mass_radius_crop


__all__ = [
    "center",
    "is_centered",
    "Centralizer",
    "potential",
    "Potentializer",
    "star_align",
    "is_star_aligned",
    "Aligner",
    "center_and_align",
    "Cutter",
    "half_star_mass_radius_crop",
]


# =============================================================================
# FUNCTIONS
# =============================================================================


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
    center = Centralizer()
    galaxy = center.transform(galaxy)  # Centering
    align = Aligner(r_cut)  # D: asi funcionaria ahpra no?
    galaxy = align.transform(galaxy)

    return galaxy


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
    # D:Creo que esto habria que cambiarlo con los checkers nuevos
    check_center = is_centered(galaxy, rtol=rtol, atol=atol)
    # check_center = pcenter.is_centered(galaxy, rtol=rtol, atol=atol)
    check_align = is_star_aligned(galaxy, r_cut=r_cut, rtol=rtol, atol=atol)

    # check_align = salign.is_star_aligned(
    #    galaxy, r_cut=r_cut, rtol=rtol, atol=atol
    # )
    return check_center and check_align

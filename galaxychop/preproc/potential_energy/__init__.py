# =============================================================================
# DOCS
# =============================================================================

"""Module for calculus of potential energy."""

# =============================================================================
# IMPORTS
# =============================================================================

import astropy.units as u
 
from attr import validators

import numpy as np

from .grispy_potential import (
    make_grid,
    potential_grispy,
)

from .. import (
    constants as const,
    core,
)

from ._base.py import (
    hparam,  # D: aca importara algun ABC de la que hereda Potential
)

try:
    from .fortran import potential as potential_f
except ImportError:  # pragma: no cover
    potential_f = None

#: The default potential backend to use.
DEFAULT_POTENTIAL_BACKEND = "numpy" if potential_f is None else "fortran"

# =============================================================================
# BACKENDS
# =============================================================================


def fortran_potential(x, y, z, m, softening):
    """
    Wrap the Fortran implementation of the gravitational potential.

    Parameters
    ----------
    x, y, z : np.ndarray
        Positions of particles. Shape: (n,1).
    m : np.ndarray
        Masses of particles. Shape: (n,1).
    softening : float, optional
        Softening parameter. Shape: (1,).

    Returns
    -------
    np.ndarray : float
        Specific potential energy of particles.

    """
    soft = np.asarray(softening)
    epot = potential_f.fortran_potential(x, y, z, m, soft)

    return epot * const.G, np.asarray


def grispy_potential(x, y, z, m, softening):
    """
    GriSPy implementation of the gravitational potential energy calculation.

    Parameters
    ----------
    x, y, z : np.ndarray
        Positions of particles. Shape: (n,1).
    m : np.ndarray
        Masses of particles. Shape: (n,1).
    softening : float, optional
        Softening parameter. Shape: (1,).

    Returns
    -------
    np.ndarray : float
        Specific potential energy of particles.

    """
    # Make the grid of the system
    l_box, grid = make_grid(x, y, z)

    # For each particle, compute its potential energy
    # Bruno:
    # Juan no me matés pls. Rev como hacer más lindo el loop.
    epot = np.empty(len(m))
    for idx, particle in enumerate(m):
        centre = np.array([x[idx], y[idx], z[idx]])
        epot[idx] = potential_grispy(
            centre,
            x,
            y,
            z,
            m,
            softening,
            bubble_size=5 * softening,
            shell_width=0.1 * l_box,
            l_box=l_box,
            grid=grid,
        )
    # Bruno:
    # Está chanchísimo escrito esto. WIP; btw, le saqué el uso de \
    # "G" adentro de la func y lo aplico desde acá directamente, \
    # apra unificar...

    return epot * const.G, np.asarray


def numpy_potential(x, y, z, m, softening):
    """
    Numpy implementation for the gravitational potential energy calculation.

    Parameters
    ----------
    x, y, z : np.ndarray
        Positions of particles. Shape: (n,1).
    m : np.ndarray
        Masses of particles. Shape:(n,1).
    softening : float, optional
        Softening parameter. Shape: (1,).

    Returns
    -------
    np.ndarray : float
        Specific potential energy of particles.

    """
    dist = np.sqrt(
        np.square(x - x.reshape(-1, 1))
        + np.square(y - y.reshape(-1, 1))
        + np.square(z - z.reshape(-1, 1))
        + np.square(softening)
    )

    np.fill_diagonal(dist, 0.0)

    flt = dist != 0
    mdist = np.divide(m, dist, out=np.zeros_like(dist), where=flt)

    return mdist.sum(axis=1) * const.G, np.asarray


# =============================================================================
# API
# =============================================================================

POTENTIAL_BACKENDS = {
    "fortran": fortran_potential,
    "grispy": grispy_potential,
    "numpy": numpy_potential,
}


class Potential:  # D: Hereda de un tentativo GalaxyTransformerABC ?
    """
    Potential energy calculation.

    Given the positions and masses of particles, calculate
    their specific gravitational potential energy.

    Parameters
    ----------
    galaxy : ``Galaxy class`` object

    backends:

    Returns
    -------
    galaxy: new ``Galaxy class`` object
        A new galaxy object with the specific potential energy of particles
        calculated.

    """

    # Idea de implementacion
    # pot = Potential(backend="frotran")
    # gal = pot.transform(gal)
    # D: no se si esto deberia ser un hiperparametro de la clase ?
    backend = hparam(
        default=POTENTIAL_BACKENDS, validator=attr.validators.instance_of(dict)
    )

    def transform(self, galaxy, *, backend):

        if galaxy.has_potential_:
            raise ValueError("galaxy potential is already calculated")

        # extract the implementation
        backend_function = POTENTIAL_BACKENDS[backend]

        # convert the galaxy in multiple arrays
        df = galaxy.to_dataframe(attributes=["x", "y", "z", "m", "softening"])
        x = df.x.to_numpy(dtype=np.float32)
        y = df.y.to_numpy(dtype=np.float32)
        z = df.z.to_numpy(dtype=np.float32)
        m = df.m.to_numpy(dtype=np.float32)
        softening = np.asarray(df.softening.max(), dtype=np.float32)

        # cleanup df
        del df

        # execute the function and return
        pot, postproc = backend_function(x, y, z, m, softening)

        # cleanup again
        del x, y, z, m, softening

        # apply the post process to the final potential
        pot = postproc(pot)

        # recreate a new galaxy
        num_s = len(galaxy.stars)
        num = len(galaxy.stars) + len(galaxy.dark_matter)

        pot_s = pot[:num_s]
        pot_dm = pot[num_s:num]
        pot_g = pot[num:]

        new = galaxy.disassemble()

        new.update(
            potential_s=-pot_s * (u.km / u.s) ** 2,
            potential_dm=-pot_dm * (u.km / u.s) ** 2,
            potential_g=-pot_g * (u.km / u.s) ** 2,
        )

        return core.mkgalaxy(**new)
# This file is part of
# the galaxy-chop project (https://github.com/vcristiani/galaxy-chop)
# Copyright (c) Cristiani, et al. 2021, 2022, 2023
# License: MIT
# Full Text: https://github.com/vcristiani/galaxy-chop/blob/master/LICENSE.txt

# =============================================================================
# DOCS
# =============================================================================

"""Utilities for align the galaxies."""

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
# Bruno: ¿?


def _make_mask(x, y, z, r_cut):
    # Bruno:
    # ¿No debería ser ~"aux_r = self.r_cut"?
    r = np.sqrt(x**2 + y**2 + z**2)

    if r_cut is None:
        return np.repeat(True, len(r))

    return np.where(r < r_cut)


def _get_rot_matrix(m, x, y, z, Jx, Jy, Jz, r_cut):
    """
    Rotation matrix calculation.

    Calculates the rotation matrix that aligns the
    TOTAL angular momentum of
    the particles with the z-axis. Optionally,
    only particles within a cutting
    radius `(r_cut)` can be used.

    Parameters
    ----------
    m : np.ndarray
        Masses of particles. Shape: (n,1).
    x, y, z : np.ndarray
        Positions x, y, z of particles. Shape: (n,1).
    Jx, Jy, Jz : np.ndarray
        Components of angular momentum of particles.
        Shape: (n,1).
    r_cut : float, optional
        Default value =  None. If it's provided,
        it must be positive and the
        rotation matrix `A` is calculated from
        the particles with radii smaller
        than r_cut.

    Returns
    -------
    A : np.ndarray
        Rotation matrix. Shape: (3,3).

    """
    mask = _make_mask(x, y, z, r_cut)

    mjx, mjy, mjz = m * Jx, m * Jy, m * Jz

    rjx = np.sum(mjx[mask])
    rjy = np.sum(mjy[mask])
    rjz = np.sum(mjz[mask])

    rjp = np.sqrt(rjx**2 + rjy**2)
    rj = np.sqrt(rjx**2 + rjy**2 + rjz**2)

    e1x = rjy / rjp
    e1y = -rjx / rjp
    e1z = 0.0

    e2x = rjx * rjz / (rjp * rj)
    e2y = rjy * rjz / (rjp * rj)
    e2z = -rjp / rj

    e3x = rjx / rj
    e3y = rjy / rj
    e3z = rjz / rj

    A = np.array(([e1x, e1y, e1z], [e2x, e2y, e2z], [e3x, e3y, e3z]))

    return A


# =============================================================================
# ALIGNER CLASS
# =============================================================================


class Aligner(GalaxyTransformerABC):
    """
    Aligner class.

    Given the positions and velocities of particles, check and
    align them to make the Z-axis perpendicular to the galaxy-plane
    i.e. make the "face-on" projection the new XY plane.

    """

    # Bruno:
    # Cambio esto porque no le gustó...
    r_cut = hparam(default=30)
    #def __init__(self, r_cut=hparam(default=30)):
    #    self.r_cut = r_cut

    # r_cut = hparam(default=30)
    # D:
    # Bruno: Nos suelen gustar ~30 kpc. Pero si 3 r_half es mucho menor,
    # hay que tener cuidado...

    @doc_inherit(GalaxyTransformerABC.transform)
    def transform(self, galaxy):
        # D: el transform no deberia ir adentro del argumento de la clase?
        # Bruno:
        # Falta acomodar/checkear que las velocidades estén corregidas por
        # v_CM para volver a calcular el Jx, Jy y Jz como corresponde (!)
        return star_align(galaxy, r_cut=self.r_cut)  # Suponiendo que así queremos...

    @doc_inherit(GalaxyTransformerABC.checker)
    def checker(self, galaxy, **kwargs):
        # Bruno:
        # ¿Mejor forma que no sea con el "**"? ¿Aclarar? Ojo...
        return is_star_aligned(galaxy, **kwargs)


# =============================================================================
# API FUNCTIONS
# =============================================================================
# Bruno: Rev como dejar este título de sec como en otros módulos...


def star_align(galaxy, *, r_cut=None):
    """
    Align the galaxy.

    Rotates the positions, velocities and angular momentum of the
    particles so that the total angular moment of the stars particles coincides
    with the z-axis. Optionally, only stars particles within a cutting radius
    `(r_cut)` can be used to calculate the rotation matrix.

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
        A new galaxy object with their total angular momentum aligned with the
        z-axis.

    """
    if r_cut is not None and r_cut <= 0.0:
        raise ValueError("r_cut must not be lower than 0.")

    # declare all the different groups of columns
    pos_columns = ["x", "y", "z"]
    vel_columns = ["vx", "vy", "vz"]

    # Now we extract only the needed column to rotate the galaxy
    # Note: for stars we need more columns to calculate the rotation matrix
    stars_df = galaxy.stars.to_dataframe(
        attributes=["m", "Jx", "Jy", "Jz"] + pos_columns + vel_columns
    )
    dm_df = galaxy.dark_matter.to_dataframe(
        attributes=pos_columns + vel_columns
    )
    gas_df = galaxy.gas.to_dataframe(attributes=pos_columns + vel_columns)

    # now we can calculate the rotation matrix
    A = _get_rot_matrix(
        m=stars_df["m"].values,
        x=stars_df["x"].values,
        y=stars_df["y"].values,
        z=stars_df["z"].values,
        Jx=stars_df["Jx"].values,
        Jy=stars_df["Jy"].values,
        Jz=stars_df["Jz"].values,
        r_cut=r_cut,
    )

    # we rotate  independently positions and velocities in stars dm and gas
    pos_rot_s = np.dot(A, stars_df[pos_columns].T.values)
    vel_rot_s = np.dot(A, stars_df[vel_columns].T.values)

    pos_rot_dm = np.dot(A, dm_df[pos_columns].T.values)
    vel_rot_dm = np.dot(A, dm_df[vel_columns].T.values)

    pos_rot_g = np.dot(A, gas_df[pos_columns].T.values)
    vel_rot_g = np.dot(A, gas_df[vel_columns].T.values)

    # recreate the valaxy
    new = galaxy.disassemble()

    new.update(
        x_s=pos_rot_s.T[:, 0],
        y_s=pos_rot_s.T[:, 1],
        z_s=pos_rot_s.T[:, 2],
        vx_s=vel_rot_s.T[:, 0],
        vy_s=vel_rot_s.T[:, 1],
        vz_s=vel_rot_s.T[:, 2],
        x_dm=pos_rot_dm.T[:, 0],
        y_dm=pos_rot_dm.T[:, 1],
        z_dm=pos_rot_dm.T[:, 2],
        vx_dm=vel_rot_dm.T[:, 0],
        vy_dm=vel_rot_dm.T[:, 1],
        vz_dm=vel_rot_dm.T[:, 2],
        x_g=pos_rot_g.T[:, 0],
        y_g=pos_rot_g.T[:, 1],
        z_g=pos_rot_g.T[:, 2],
        vx_g=vel_rot_g.T[:, 0],
        vy_g=vel_rot_g.T[:, 1],
        vz_g=vel_rot_g.T[:, 2],
    )

    return data.mkgalaxy(**new)


def is_star_aligned(galaxy, *, r_cut=None, rtol=1e-05, atol=1e-08):
    """
    Validate if the galaxy is aligned.

    Parameters
    ----------
    galaxy : ``Galaxy class`` object
    r_cut : float, optional
        Default value =  None. If it's provided,
        it must be positive and the
        rotation matrix `A` is calculated from
        the particles with radii smaller
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
    # Now we extract only the needed column to rotate the galaxy
    df = galaxy.stars.to_dataframe(
        attributes=["m", "x", "y", "z", "Jx", "Jy", "Jz"]
    )

    mask = _make_mask(df.x.values, df.y.values, df.z.values, r_cut)

    Jxtot = np.sum(df.Jx.values[mask] * df.m.values[mask])
    Jytot = np.sum(df.Jy.values[mask] * df.m.values[mask])
    Jztot = np.sum(df.Jz.values[mask] * df.m.values[mask])
    Jtot = np.sqrt(Jxtot**2 + Jytot**2 + Jztot**2)

    return np.allclose(Jztot, Jtot, rtol=rtol, atol=atol)

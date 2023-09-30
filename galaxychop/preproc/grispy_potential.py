#Bruno:
# Para comenzar a hacer algo, veamos la implementación \
# del potencial con GriSPy (aunqeu ya esté testeado que \
# no es el algoirtmo más rápido ni más preciso...).

#Bruno:
# OJO que el copy de esto es sobre la branch "master" \
# =/= "dev", así que puede que algo esté desfazado...

#Bruno:
# Estaría bueno tener una clase que abarque los métodos \
# para calcular potencial, así como los descomponedores \
# Gausseanos, etc...; Por ahora, me limito a copiar la \
# metodología de las funciones anteriores de cálculo de \
# potencial ('numpy' y 'fortran')...

# This file is part of
# the galaxy-chop project (https://github.com/vcristiani/galaxy-chop)
# Copyright (c) 2021, Valeria Cristiani
# License: MIT
# Full Text: https://github.com/vcristiani/galaxy-chop/blob/master/LICENSE.txt

# =============================================================================
# DOCS
# =============================================================================

"""Utilities to make a grid and run Nearest Neighbor Searchs via GriSPy."""

# =============================================================================
# IMPORTS
# =============================================================================

import grispy as gsp

import numpy as np

# =============================================================================
# BACKENDS
# =============================================================================


def make_grid(x, y, z, n_cells=2**4):
    """
    Space grid maker.

    Make a grid of the volume and index galaxy's particles to each cell.

    Parameters
    ----------
    x, y, z : np.ndarray
        Positions of particles. Shape: (n,1).
    n_cells : float
        Number of cells that makes the grid. Shape: (1,).
        #Bruno:
        # Aclarar el valor default*

    Returns
    ------- *Cambiar
    l_box : float
        Length of the size of the smalles cube
          that encloses all particles. Shape: (1,).
    #Bruno:
    # ¿Cómo se pone bien? Revisar...
    grid : ``GriSPy`` object
        Grid populated with all the galaxy particles.
          Shape: (n,1).
    #Bruno:
    # ¿Por qué me corre 2 espacios en docs? Revsiar...

    """
    # Size of the box that contains all particles
    l_box = max(np.abs([max(x)-min(x),
                       max(y)-min(y),
                       max(z)-min(z)]))

    # Make the grid (n_cells ~ 2**4 works well for 1e+4 ~ 1e+5 particles)
    grid = gsp.GriSPy(x, y, z, n_cells)

    #Bruno:
    # No uso la masa "m", pero ojo con los procedures en relación a los
    # demás métodos...

    return l_box, grid


def potential_grispy(centre, x, y, z, m, softening,
                     bubble_size, shell_width, l_box, grid):
    """
    Compute the potential of a particle given the grid and the particles.

    Given the particle to compute its potential energy, iteratively
    make shells to aproximate their monopole contribution to the total
    potential energy.
    #Bruno:
    # Pero claro, todas las otras funciones vienen comiéndose "galaxias",
    ¿Debería mantenerme con eso acá? Estas funciones van de la mano,
    así que ojo con mi implementación naïve...

    Parameters
    ---------- *Revisar
    centre : np.array
        3D spatial position of the particle to compute its potential.
        Shape: (n,3).
    x,y,z : np.array
        3D spatial position of all the galaxy particles. Shape: (n,3).
    m : np.array
        Individual masses of all the galaxy particles. Shape: (n,1).
    bubble_size : float
        Radii of the sphere that will contain the closest particles,
        which potential contribution will be calculated via direct-
        sumation. Shape: (1,).
    shell_width : float
        Width of the consecutive shells that will contain further particles.
        The shell's potential contribution will be aproximated by its
        monopole term. Shape: (1,).
    l_box : float
        Size of the box that contains all particles. This defines the upper
        limit value of the shells to implement the GriSPy's NNS. Shape: (1,).
    grid : ``GriSPy`` object
        Spatial grid populated by the galaxy particles, needed to do the NNS.

    Returns
    -------
    pot_shells : float
        Potential of the given particle through the shells' monopole
        aproximation. Shape: (1,).

    """
    #Bruno:
    # En este caso no hace falta el potencial ¿O sí?
    # Lo cambié para que, si ya tiene potencial, largue error
    #if galaxy.has_potential_:
    #    raise ValueError("galaxy already has the potential energy")
    #Bruno:
    # Entonces sí debería comer galaxias...

    # Use the bubble method to find the closest particles
    bubble_dist, bubble_ind = grid.bubble_neighbors(
        centre,
        distance_upper_bound=bubble_size
    )

    # Compute the potential contribution via direct-sumation
    # of these bubble's particles
    pot_shells = 0.  # The potential variable.
    for idx, distance in enumerate(bubble_dist[0]):
        if distance > 0.:
            pot_shells -= m[bubble_ind[0][idx]]/distance
        else:
            continue

        #Bruno:
        # Otra versión de esto (mucho más lenta, pero creo que correcta)
        #d_and_soft = np.sqrt(np.square(distance) + np.square(softening))
        #pot_shells -= m[bubble_ind[0][idx]]/d_and_soft

    d_min_shell = bubble_size  # Shell's lower limit to initialize the loop
    while d_min_shell < l_box:
        # Use the shell method to populate it
        shell_dist, shell_ind = grid.shell_neighbors(
            centre,
            distance_lower_bound=d_min_shell,
            distance_upper_bound=d_min_shell + shell_width
        )

        # Compute the monopole potential contribution of this shell
        for idx, distance in enumerate(shell_dist[0]):
            # Due to non-periodicity, distance > 0 always
            pot_shells -= m[shell_ind[0][idx]]/distance

            # Otra versión de esto (mucho más lenta, pero creo que correcta)
            #d_and_soft = np.sqrt(np.square(distance) + np.square(softening))
            #pot_shells -= m[shell_ind[0][idx]]/d_and_soft

        d_min_shell += shell_width  # Repeat for next shell (further away)

    return pot_shells

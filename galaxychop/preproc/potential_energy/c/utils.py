# This file is part of
# the galaxy-chop project (https://github.com/vcristiani/galaxy-chop)
# Copyright (c) Cristiani, et al. 2021, 2022, 2023
# License: MIT
# Full Text: https://github.com/vcristiani/galaxy-chop/blob/master/LICENSE.txt

# Bruno:
# Este archivo vino junto al potential.c, así que no me animé \
# a tocar nada...; Falta integrarlo, ver cómo llamarlo en el \
# "pot_nrg.py" y testearlo, claro. ¿Hace falta que lo "haga \
# lindo", o es mejor rehacerlo todo?

# Bruno:
# Si bien lo voy a dejar en esta carpeta de momento, lo mejor \
# es llevarlo a otra carpeta "potentials" junto con lo de \
# fortran y lo de GriSPy...
"""C implementations."""

import os
from ctypes import c_int

import numpy as np
import numpy.ctypeslib as npct


def check(exe, recompilar=True):
    """Check function."""
    if not os.path.exists("%s.so" % exe) or recompilar:
        print("Compile %s.c" % exe)
        comando = "gcc -c -O3 -fPIC -Wall -lm -fopenmp %s.c -o %s.o" % (
            exe,
            exe,
        )
        os.system(comando)
        if os.path.exists("%s.o" % exe):
            print("Compile %s.so" % exe)
            comando = "gcc -shared -Wall -lm -fopenmp %s.o -o %s.so" % (
                exe,
                exe,
            )
            os.system(comando)
        else:
            print("Error!!!\n")
            exit(0)

    return


def calcula_potencial(Npart, mp, x, y, z):
    """Calcula potential function."""
    array_1d_float = npct.ndpointer(
        dtype=np.float32, ndim=1, flags="C_CONTIGUOUS"
    )
    # D: no se usa
    # array_1d_int = npct.ndpointer(dtype=np.int32, ndim=1,
    #  flags="C_CONTIGUOUS")

    libcd = npct.load_library("potential.so", ".")

    libcd.calculate_potential.restype = None
    libcd.calculate_potential.argtypes = [
        c_int,
        array_1d_float,
        array_1d_float,
        array_1d_float,
        array_1d_float,
        array_1d_float,
    ]

    Ep = np.zeros(Npart, dtype=np.float32)

    libcd.calculate_potential(Npart, mp, x, y, z, Ep)

    return Ep

# Revisando las slides de Juan, encontré que hay otra salida
# a "traducir" el Octree de C a Cython -> ¡Un wrapper!
# Veamos qué podemos/debemos hacer...

cdef extern from "potential.c":
    void calculate_potential(
        const int npart, const float *mp, const float *x, const float *y, const float *z, float *Ep
        )

# Ojo con los pointers, ¿molestarán a la hora de usar el Potentializer?
# En https://cython.readthedocs.io/en/latest/src/userguide/external_C_code.html
# habla de usar "cdef ..." y "&" para meter estas variables... ¿?

def py_calculate_potential():
    return calculate_potential(npart,mp,x,y,z,Ep)


# Dejemos esto. Probamos hacer un wrap de Python en C (new folder)
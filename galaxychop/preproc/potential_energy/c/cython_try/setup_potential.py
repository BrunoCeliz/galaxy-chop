# Siguiendo el walkthrough para realizar el wrapp,
# copio lo de
# https://stavshamir.github.io/python/making-your-c-library-callable-from-python-by-wrapping-it-with-cython/

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

examples_extension = Extension(
    name="octree_pot",
    sources=["potential_cwrapperpyexamples.pyx"],
    libraries=["potential"],
    library_dirs=["lib"],
    include_dirs=["lib"],
)
setup(name="octree_pot", ext_modules=cythonize([examples_extension]))

# To do:
# python setup_potential.py build_ext --inplace
# ¿Hacer cosas dentro de un "makefile"?

# Para llamarla:
# import pyexamples
# if __name__ == '__main__':
#     octree_pot.py_calculate_potential(args)

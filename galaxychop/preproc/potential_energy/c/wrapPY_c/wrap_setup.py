from distutils.core import setup, Extension

setup(
    name="fputs",
    ext_modules=[Extension("potential_c_module", ["potential.c"])],
)

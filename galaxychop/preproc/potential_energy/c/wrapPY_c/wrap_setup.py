# from distutils.core import setup
# from distutils.extension import Extension
# from Pyrex.Distutils import build_ext

from setuptools import setup, find_packages, Extension

setup(
    name='pot',
    version='0.1.0',
    packages=find_packages(),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Natural Language :: English',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: Implementation :: CPython',
    ],
    ext_modules=[
        Extension(
            # the qualified name of the extension module to build
            'pot.pot',
            # the files to compile into our module relative to ``setup.py``
            ['potential_wrapped.c'],
        ),
    ],
)
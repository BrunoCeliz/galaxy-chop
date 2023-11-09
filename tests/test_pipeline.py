# This file is part of
# the galaxy-chop project (https://github.com/vcristiani/galaxy-chop)
# Copyright (c) Cristiani, et al. 2021, 2022, 2023
# License: MIT
# Full Text: https://github.com/vcristiani/galaxy-chop/blob/master/LICENSE.txt


# =============================================================================
# DOCS
# =============================================================================

"""test for galaxy-chop.pipelines

"""


# =============================================================================
# IMPORTS
# =============================================================================


# import galaxychop as gchop

from galaxychop import pipeline
from galaxychop.models.gaussian_mixture import GaussianMixture
from galaxychop.models.threshold import JThreshold
from galaxychop.preproc.pcenter import Centralizer
from galaxychop.preproc.potential_energy import Potentializer
from galaxychop.preproc.salign import Aligner

import pytest

# =============================================================================
# TESTS
# =============================================================================

# D: entiendo que esto checkea que tan bien esta
# funcionando el mkpipe utilizando los metodos
# no el pipeline


def test_pipeline_mkpipe(read_hdf5_galaxy):
    gal = read_hdf5_galaxy("gal394242.h5")

    steps = [
        Centralizer(),
        Aligner(),
        JThreshold(),
    ]

    expected = gal
    # D:
    # center = gchop.preproc.Centralizer()
    # gal = center.transform(gal_in)
    # align = gchop.preproc.Aligner()
    # expected_transform = align.transform(gal)

    # decomp = gchop.models.JThreshold()
    # expected_components = decomp.decompose(expected_transform)

    # D: lo anterior lo hice a mano pero entiendo que
    # hace lo mismo que esto
    for step in steps[:-1]:
        expected = step.transform(expected)
    expected = steps[-1].decompose(expected)

    pipe = pipeline.mkpipe(*steps)
    result = pipe.decompose(gal)

    assert result.values_equals(expected)
    assert len(pipe) == len(steps)
    assert steps == [s for _, s in pipe.steps]
    for s in pipe.named_steps.values():
        assert s in steps


def test_pipeline_slicing():
    steps = [
        Centralizer(),
        Aligner(),
        Potentializer(),
        JThreshold(),
    ]

    pipe = pipeline.mkpipe(*steps)

    for idx, step in enumerate(steps):
        assert pipe[idx] == step

    for name, step in pipe.named_steps.items():
        assert pipe[name] == step

    assert [s for _, s in pipe[2:].steps] == steps[2:]

    with pytest.raises(ValueError):
        pipe[::2]

    with pytest.raises(KeyError):
        pipe[None]


def test_pipeline_not_transformer_fail():
    steps = [JThreshold(), GaussianMixture()]
    with pytest.raises(TypeError):
        pipeline.mkpipe(*steps)


def test_pipeline_not_decomposer_fail():
    steps = [Centralizer(), Aligner()]
    with pytest.raises(TypeError):
        pipeline.mkpipe(*steps)


def test_pipeline_name_not_str():
    with pytest.raises(TypeError):
        pipeline.GchopPipeline(
            steps=[(..., Centralizer()), ("final", JThreshold())]
        )
    with pytest.raises(TypeError):
        pipeline.GchopPipeline(
            steps=[("first", Centralizer()), (..., JThreshold())]
        )


# D:
# Me queda pendiente poder ver una forma de comparar resultados
# tengo que ver bien como lo hacen los otros test
# para poder comparar lo que va haciendo los metodos
# o simplemente confio en los test de los metodos

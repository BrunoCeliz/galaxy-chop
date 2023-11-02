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

import pytest

from galaxychop import pipeline
from galaxychop.preproc.pcenter import Centralizer
from galaxychop.preproc.saling import Aligner
from galaxychop.preproc.potential_energy import Potential
from galaxychop.models.gaussian_mixture import GaussianMixture 
# from galaxychop.models.histogram import JEHistogram, JHistogram
# from galaxychop.models.kmeans import KMeans
from galaxychop.models.threshold import JThreshold

# =============================================================================
# TESTS
# =============================================================================

# D: tengo duda de aca como comparar el assert
# def test_pipeline_mkpipe(decision_matrix):
#    dm = decision_matrix(seed=42)

#    steps = [
#        g(),
#        StandarScaler(target="matrix"),
#        CRITIC(correlation="spearman"),
#        CRITIC(),
#        TOPSIS(),
#    ]

#    expected = dm
#    for step in steps[:-1]:
#        expected = step.transform(expected)
#    expected = steps[-1].evaluate(expected)

#    pipe = pipeline.mkpipe(*steps)
#    result = pipe.evaluate(dm)

#    assert result.values_equals(expected)
#    assert len(pipe) == len(steps)
#    assert steps == [s for _, s in pipe.steps]
#    for s in pipe.named_steps.values():
#        assert s in steps

# D: estos trcuatro es creo que deberan andar

def test_pipeline_slicing():
    steps = [
        Centralizer(),
        Aligner(),
        Potential(),
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
    steps = [Centralizer(),Aligner()]
    with pytest.raises(TypeError):
        pipeline.mkpipe(*steps)


def test_pipeline_name_not_str():
    with pytest.raises(TypeError):
        pipeline.GchopPipeline(steps=[(..., Centralizer()), ("final", JThreshold())])
    with pytest.raises(TypeError):
        pipeline.GchopPipeline(steps=[("first", Centralizer()), (..., JThreshold())])
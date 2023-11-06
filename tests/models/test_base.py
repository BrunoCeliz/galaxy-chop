# This file is part of
# the galaxy-chop project (https://github.com/vcristiani/galaxy-chop)
# Copyright (c) Cristiani, et al. 2021, 2022, 2023
# License: MIT
# Full Text: https://github.com/vcristiani/galaxy-chop/blob/master/LICENSE.txt


import galaxychop as gchop

import numpy as np

import pandas as pd

import pytest


# =============================================================================
# COMPONENTS
# =============================================================================


@pytest.mark.model
@pytest.mark.parametrize("probs", [True, False])
def test_Components(probs):
    random = np.random.default_rng(42)

    labels = random.integers(0, 3, 100)
    ptypes = np.ones(100)
    probabilities = random.normal(size=100) if probs else None
    mass = random.normal(size=100)

    components = gchop.models.Components(
        labels=labels,
        ptypes=ptypes,
        probabilities=probabilities,
        m=mass,
        lmap={},
    )

    assert len(components) == 100

    expected_repr = (
        "<Components length=100, labels={0, 1, 2}, "
        f"probabilities={probs}, lmap=False>"
    )
    assert repr(components) == expected_repr


@pytest.mark.model
@pytest.mark.parametrize("probs", [True, False])
def test_Components_bad_len(probs):
    random = np.random.default_rng(42)

    labels = random.integers(0, 3, 100)
    ptypes = np.ones(99)
    mass = random.normal(size=95)
    probabilities = random.normal(size=98) if probs else None

    with pytest.raises(ValueError):
        gchop.models.Components(
            labels=labels,
            ptypes=ptypes,
            probabilities=probabilities,
            m=mass,
            lmap={},
        )


@pytest.mark.model
@pytest.mark.parametrize("probs", [True, False])
def test_Components_to_dataframe(probs):
    random = np.random.default_rng(42)

    labels = random.integers(0, 3, 100)
    ptypes = np.ones(100)
    mass = random.normal(size=100)
    probabilities = random.normal(size=100) if probs else None

    components = gchop.models.Components(
        labels=labels,
        ptypes=ptypes,
        probabilities=probabilities,
        m=mass,
        lmap={},
    )

    expected = pd.DataFrame(
        {
            "m": mass,
            "labels": labels,
            "ptypes": ptypes,
            "lmap": labels.astype(object),
        }
    )

    if probs:
        probs_df = pd.DataFrame({"probs_0": probabilities})
        expected = pd.concat([expected, probs_df], axis=1)

    pd.testing.assert_frame_equal(components.to_dataframe(), expected)


@pytest.mark.model
@pytest.mark.parametrize("probs", [True, False])
def test_Components_describe(probs):
    random = np.random.default_rng(42)

    labels = random.integers(0, 3, 100)
    ptypes = np.ones(100)
    mass = random.normal(loc=1000933.2, scale=252304.96, size=100)
    probabilities = random.uniform(size=(100, 3)) if probs else None

    components = gchop.models.Components(
        labels=labels,
        ptypes=ptypes,
        probabilities=probabilities,
        m=mass,
        lmap={},
    )

    expected_dict = {
        ("Particles", "Size"): {0: 27, 1: 33, 2: 40},
        ("Particles", "Fraction"): {0: 0.27, 1: 0.33, 2: 0.4},
        ("Deterministic mass", "Size"): {
            0: 26087288.44203574,
            1: 32763518.959382985,
            2: 38057037.01262353,
        },
        ("Deterministic mass", "Fraction"): {
            0: 0.26919687048838753,
            1: 0.33808944113336864,
            2: 0.39271368837824394,
        },
    }

    if probs:
        expected_dict.update(
            {
                ("Probabilistic mass", "Size"): {
                    0: 50689931.17716688,
                    1: 46151058.83650705,
                    2: 48717115.96184553,
                },
                ("Probabilistic mass", "Fraction"): {
                    0: 0.523073560078298,
                    1: 0.4762365638773781,
                    2: 0.5027159179570635,
                },
            }
        )

    expected = pd.DataFrame.from_dict(expected_dict)

    result = components.describe()

    pd.testing.assert_frame_equal(result, expected)


# =============================================================================
# DECOMPOSER ABC
# =============================================================================


@pytest.mark.model
def test_GalaxyDecomposerABC_not_implemethed():
    class Decomposer(gchop.models.GalaxyDecomposerABC):
        def get_attributes(self):
            return super().get_attributes()

        def split(self, X, y, attributes):
            return super().split(X, y, attributes)

        def get_rows_mask(self, X, y, attributes):
            return super().get_rows_mask(X, y, attributes)

    decomposer = Decomposer()

    with pytest.raises(NotImplementedError):
        decomposer.get_attributes()

    with pytest.raises(NotImplementedError):
        decomposer.split(None, None, None)

    with pytest.raises(NotImplementedError):
        decomposer.get_rows_mask(None, None, None)


@pytest.mark.model
@pytest.mark.parametrize(
    "bins_value", [None, (1.0,), (1.0, 2.0, 3.0), (1.0, 2)]
)
def test_GalaxyDecomposerABC_invalid_bins(bins_value):
    class Decomposer(gchop.models.GalaxyDecomposerABC):
        def get_attributes(self):
            ...

        def split(self, X, y, attributes):
            ...

        def get_rows_mask(self, X, y, attributes):
            ...

    with pytest.raises(ValueError):
        Decomposer(cbins=bins_value)


@pytest.mark.model
def test_GalaxyDecomposerABC_repr():
    class Decomposer(gchop.models.GalaxyDecomposerABC):
        other = gchop.models.hparam(default=1)

        def get_attributes(self):
            return ["normalized_star_energy", "eps", "eps_r"]

        def split(self, X, y, attributes):
            ...

        def get_rows_mask(self, X, y, attributes):
            ...

    decomposer = Decomposer(cbins=(0.3, 0.2), reassign=True, other="zaraza")
    result = repr(decomposer)
    expected = "Decomposer(cbins=(0.3, 0.2), reassign=True, other='zaraza')"

    assert result == expected


@pytest.mark.model
def test_GalaxyDecomposerABC_attributes_matrix(read_hdf5_galaxy):
    gal = read_hdf5_galaxy("gal394242.h5")
    gal = gchop.preproc.star_align(gchop.preproc.center(gal))

    class Decomposer(gchop.models.GalaxyDecomposerABC):
        def get_attributes(self):
            ...

        def split(self, X, y, attributes):
            ...

        def get_rows_mask(self, X, y, attributes):
            ...

    decomposer = Decomposer()

    attributes = ["x", "eps"]

    X, t = decomposer.attributes_matrix(gal, attributes=attributes)

    # check types stars-dm-gas
    assert np.all(t[: len(gal.stars)] == gchop.ParticleSetType.STARS.value)
    assert np.all(
        t[len(gal.stars) : len(gal.dark_matter)]  # noqa
        == gchop.ParticleSetType.DARK_MATTER.value
    )
    assert np.all(
        t[len(gal.stars) + len(gal.dark_matter) : len(gal.gas)]  # noqa
        == gchop.ParticleSetType.GAS.value
    )

    # check as_dataframe attrs
    assert np.all(X[:, 0] == gal.to_dataframe(attributes=["x"])["x"])

    # check jcirc eps
    jcirc = gal.stellar_dynamics()

    X_stars = X[t == gchop.ParticleSetType.STARS.value]
    assert np.array_equal(X_stars[:, 1], jcirc.eps, equal_nan=True)

    X_nostars = X[t != gchop.ParticleSetType.STARS.value]
    assert np.all(np.isnan(X_nostars[:, 1]))


@pytest.mark.model
def test_GalaxyDecomposerABC_complete_labels():
    class Decomposer(gchop.models.GalaxyDecomposerABC):
        def get_attributes(self):
            ...

        def split(self):
            ...

        def get_rows_mask(self, X, y, attributes):
            ...

    decomposer = Decomposer()

    X = np.random.rand(3, 4)
    labels = [1, 1]
    rows_mask = [True, False, True]

    result = decomposer.complete_labels(
        X=X, labels=labels, rows_mask=rows_mask
    )

    assert np.array_equal(result, [1, np.nan, 1], equal_nan=True)


@pytest.mark.model
def test_GalaxyDecomposerABC_decompose(read_hdf5_galaxy):
    gal = read_hdf5_galaxy("gal394242.h5")
    gal = gchop.preproc.star_align(gchop.preproc.center(gal))

    class Decomposer(gchop.models.GalaxyDecomposerABC):
        def get_attributes(self):
            return ["x"]

        def split(self, X, y, attributes):
            return np.full(len(X), 100), None

        def get_rows_mask(self, X, y, attributes):
            return y == 2

    decomposer = Decomposer()

    components = decomposer.decompose(gal)

    assert (components.ptypes == "stars").sum() == len(gal.stars)
    assert (components.ptypes == "dark_matter").sum() == len(gal.dark_matter)
    assert (components.ptypes == "gas").sum() == len(gal.gas)

    assert np.all(components.labels[components.ptypes == "gas"] == 100)
    assert np.all(np.isnan(components.labels[components.ptypes != "gas"]))

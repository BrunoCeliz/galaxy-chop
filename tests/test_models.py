# This file is part of
# the galxy-chop project (https://github.com/vcristiani/galaxy-chop)
# Copyright (c) 2020, Valeria Cristiani
# License: MIT
# Full Text: https://github.com/vcristiani/galaxy-chop/blob/master/LICENSE.txt

"""Test dynamical decomposition methods."""

# =============================================================================
# IMPORTS
# =============================================================================

from galaxychop import models

import numpy as np

import pytest

from sklearn.cluster import KMeans

# =============================================================================
# TESTS
# =============================================================================


def test_GCKmeans(mock_real_galaxy):
    """Test GCKmeans."""
    gal = mock_real_galaxy

    gckmeans = models.GCKmeans(n_clusters=5, random_state=0)
    result = gckmeans.decompose(gal)

    kmeans = KMeans(n_clusters=5, random_state=0)
    X, y = gal.values()
    expected = kmeans.fit_transform(X, y)

    np.testing.assert_allclose(result, expected, rtol=1e-9, atol=1e-06)


@pytest.mark.parametrize(
    "type_values",
    [
        "string",
        np.random.rand(1),
        np.inf,
        np.nan,
    ],
)
def test_type_error_GCDecomposeMixin_class(type_values):
    """Test type error GCDecomposeMixin."""
    with pytest.raises(TypeError):
        models.GCDecomposeMixin(type_values)


def test_GCAbadi_len(mock_real_galaxy):
    """Test the lengths of labels."""
    gal = mock_real_galaxy
    X, y = gal.values()
    abadi = models.GCAbadi(seed=10)
    abadi.decompose(gal)

    longitude = len(abadi.labels_)
    assert np.shape(X) == (longitude, 10)


def test_GCAbadi_outputs(mock_real_galaxy):
    """Test outputs of GCAbadi model."""
    gal = mock_real_galaxy
    abadi = models.GCAbadi(seed=10)
    abadi.decompose(gal)

    labels = abadi.labels_
    (comp0,) = np.where(labels == 0)
    (comp1,) = np.where(labels == 1)
    (comp_nan,) = np.where(labels == -1)
    len_lab = len(labels[comp0]) + len(labels[comp1]) + len(labels[comp_nan])

    assert (labels >= -1).all() and (labels <= 1).all()
    assert len_lab == len(labels)


def test_GCAbadi_histogram(mock_real_galaxy):
    """Test the number of particles per bin."""
    gal = mock_real_galaxy
    X, y = gal.values()
    abadi = models.GCAbadi(seed=10)
    abadi.decompose(gal)
    labels = abadi.labels_
    (comp0,) = np.where(labels == 0)
    (comp1,) = np.where(labels == 1)

    full_histogram = np.histogram(X[:, 8], bins=100, range=(-1.0, 1.0))
    comp0_histogram = np.histogram(X[:, 8][comp0], bins=100, range=(-1.0, 1.0))
    comp1_histogram = np.histogram(X[:, 8][comp1], bins=100, range=(-1.0, 1.0))

    comp0_hist_plus_comp1_hist = comp0_histogram[0] + comp1_histogram[0]

    np.testing.assert_equal(comp0_hist_plus_comp1_hist, full_histogram[0])
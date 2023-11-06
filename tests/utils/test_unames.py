# This file is part of
# the galaxy-chop project (https://github.com/vcristiani/galaxy-chop)
# Copyright (c) Cristiani, et al. 2021, 2022, 2023
# License: MIT
# Full Text: https://github.com/vcristiani/galaxy-chop/blob/master/LICENSE.txt

<<<<<<< HEAD
# =============================================================================
# DOCS
# =============================================================================

"""test for galaxychop.utils.unames

"""

=======
# WIP
>>>>>>> 9bcc48e43ed744e71ca9f65f6947d3f1606d7490

# =============================================================================
# IMPORTS
# =============================================================================
<<<<<<< HEAD

import pytest

from galaxychop.utils import unames


# =============================================================================
# TEST CLASSES
# =============================================================================


def test_unique_names():
    names, elements = ["foo", "faa"], [0, 1]
    result = dict(unames.unique_names(names=names, elements=elements))
    expected = {"foo": 0, "faa": 1}
    assert result == expected


def test_unique_names_with_duplticates():
    names, elements = ["foo", "foo"], [0, 1]
    result = dict(unames.unique_names(names=names, elements=elements))
    expected = {"foo_1": 0, "foo_2": 1}
    assert result == expected


def test_unique_names_with_different_len():
    names, elements = ["foo", "foo"], [0]
    with pytest.raises(ValueError):
        unames.unique_names(names=names, elements=elements)
=======
>>>>>>> 9bcc48e43ed744e71ca9f65f6947d3f1606d7490

import pytest
from pyhyperiso.core.Math.matrix import *

def test_create_identity_matrix():
    indices = [0, 1, 2]
    identity = create_identity_matrix(indices)
    for i in indices:
        for j in indices:
            expected = 1.0 if i == j else 0.0
            assert get_element(identity, i, j) == expected

def test_get_set_element():
    mat = {}
    set_element(mat, 0, 1, 2.5)
    assert get_element(mat, 0, 1) == 2.5
    set_element(mat, 0, 1, 0.0)
    assert get_element(mat, 0, 1) == 0.0

def test_get_diagonal_elements():
    mat = {}
    for i in range(3):
        set_element(mat, i, i, i + 1)
    diag = get_diagonal_elements(mat)
    assert sorted(diag) == [0, 1, 2]

def test_invert_matrix():
    indices = [0, 1]
    mat = create_identity_matrix(indices)
    inv = invert_matrix(mat, indices)
    for i in indices:
        for j in indices:
            expected = 1.0 if i == j else 0.0
            assert get_element(inv, i, j) == expected

def test_print_matrix(capsys):
    indices = [0, 1]
    mat = create_identity_matrix(indices)
    print_matrix(mat, indices)
    captured = capsys.readouterr()
    lines = captured.out.strip().split("\n")
    assert lines == ["1.0 0.0", "0.0 1.0"]
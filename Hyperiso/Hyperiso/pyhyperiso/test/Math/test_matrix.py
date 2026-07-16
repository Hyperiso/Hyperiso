from pyhyperiso.core.Math.SparseMatrix import SparseMatrix


def test_create_identity_matrix():
    indices = [0, 1, 2]
    mat = SparseMatrix()
    for i in indices:
        mat.set(i, i, 1.0)

    for i in indices:
        for j in indices:
            expected = 1.0 if i == j else 0.0
            assert mat.get(i, j) == expected


def test_get_set_element():
    mat = SparseMatrix()
    mat.set(0, 1, 2.5)
    assert mat.get(0, 1) == 2.5

    mat.set(0, 1, 0.0)
    assert mat.get(0, 1) == 0.0


def test_get_diagonal_elements():
    mat = SparseMatrix()
    for i in range(3):
        mat.set(i, i, i + 1)

    diag = mat.diagonal()
    assert sorted(diag) == [0, 1, 2]


def test_invert_matrix():
    indices = [0, 1]
    mat = SparseMatrix()
    mat.set(0, 0, 1.0)
    mat.set(1, 1, 1.0)

    inv = mat.invert(indices)
    for i in indices:
        for j in indices:
            expected = 1.0 if i == j else 0.0
            assert inv.get(i, j) == expected

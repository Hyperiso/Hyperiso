from pyhyperiso.phyperiso.pyhyperiso import math as mb

# Matrix functions
matrix = mb.matrix

def get_diagonal_elements(matrix_data):
    """Extract diagonal indices from a sparse matrix.

    Args:
        matrix_data: Sparse matrix in dictionary form with (row, col) keys.

    Returns:
        List of diagonal indices.
    """
    return matrix.getDiagonalElements(matrix_data)

def create_identity_matrix(indices):
    """Create a sparse identity matrix.

    Args:
        indices: List of row/column indices.

    Returns:
        Dictionary representing the sparse identity matrix.
    """
    return matrix.createIdentityMatrix(indices)

def get_element(matrix_data, row, col):
    """Get value at a specific matrix position.

    Args:
        matrix_data: Sparse matrix.
        row: Row index.
        col: Column index.

    Returns:
        Value at (row, col), or 0.0 if not present.
    """
    return matrix.getElement(matrix_data, row, col)

def set_element(matrix_data, row, col, value):
    """Set value in a sparse matrix.

    Args:
        matrix_data: Sparse matrix.
        row: Row index.
        col: Column index.
        value: Value to set. If 0.0, the entry is removed.
    """
    matrix.setElement(matrix_data, row, col, value)

def invert_matrix(matrix_data, indices):
    """Invert a sparse matrix.

    Args:
        matrix_data: Sparse matrix.
        indices: List of indices defining matrix order.

    Returns:
        Inverted sparse matrix.
    """
    return matrix.invertMatrix(matrix_data, indices)

def print_matrix(matrix_data, indices):
    """Print a sparse matrix as a dense 2D table.

    Args:
        matrix_data: Sparse matrix.
        indices: List of row/column indices for ordering.
    """
    matrix.printMatrix(matrix_data, indices)
"""Sparse matrix wrapper used by low-level math utilities."""

from pyhyperiso.phyperiso.pyhyperiso import math as mb


class SparseMatrix:
    """Python wrapper around the C++ sparse matrix class.

    The wrapper exposes only a minimal API: element assignment/access, inversion
    on a selected index set, and inspection of diagonal entries. It is intended
    for advanced users who need direct access to sparse matrix primitives in the
    C++ math module.

    Example:
        >>> M = SparseMatrix()
        >>> M[0, 1] = 2.0
        >>> M[0, 1]
        2.0
    """

    def __init__(self):
        """Create an empty sparse matrix backed by C++."""
        self._cpp_obj = mb.matrix.SparseMatrix()

    def set(self, row, col, val):
        """Set one matrix entry.

        Args:
            row: Row index.
            col: Column index.
            val: Numeric value assigned to ``(row, col)``.
        """
        self._cpp_obj.set_element(row, col, val)

    def get(self, row, col):
        """Return one matrix entry.

        Args:
            row: Row index.
            col: Column index.

        Returns:
            float: Stored value returned by the C++ sparse matrix.
        """
        return self._cpp_obj.get_element(row, col)

    def print(self, indices):
        """Print a selected sparse submatrix using the C++ debug printer.

        Args:
            indices: Indices to print. The exact accepted container mirrors the
                C++ binding.
        """
        self._cpp_obj.print(indices)

    def invert(self, indices):
        """Invert a selected sparse submatrix.

        Args:
            indices: Row/column indices defining the submatrix to invert.

        Returns:
            SparseMatrix: Python wrapper around the inverted C++ matrix.
        """
        inv_cpp = self._cpp_obj.invert(indices)
        inv = SparseMatrix()
        inv._cpp_obj = inv_cpp
        return inv

    def diagonal(self):
        """Return the sparse matrix diagonal entries.

        Returns:
            Any: Container returned by the C++ binding.
        """
        return self._cpp_obj.get_diagonal_elements()

    def __getitem__(self, key):
        """Return ``self[row, col]``."""
        row, col = key
        return self.get(row, col)

    def __setitem__(self, key, value):
        """Set ``self[row, col] = value``."""
        row, col = key
        self.set(row, col, value)


__all__ = ["SparseMatrix"]

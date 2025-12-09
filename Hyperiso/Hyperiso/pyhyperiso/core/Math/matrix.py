from pyhyperiso.phyperiso.pyhyperiso import math as mb


class SparseMatrix:
    def __init__(self):
        self._cpp_obj = mb.matrix.SparseMatrix()

    def set(self, row, col, val):
        self._cpp_obj.set_element(row, col, val)

    def get(self, row, col):
        return self._cpp_obj.get_element(row, col)

    def print(self, indices):
        self._cpp_obj.print(indices)

    def invert(self, indices):
        inv_cpp = self._cpp_obj.invert(indices)
        inv = SparseMatrix()
        inv._cpp_obj = inv_cpp
        return inv

    def diagonal(self):
        return self._cpp_obj.get_diagonal_elements()

    def __getitem__(self, key):
        row, col = key
        return self.get(row, col)

    def __setitem__(self, key, value):
        row, col = key
        self.set(row, col, value)
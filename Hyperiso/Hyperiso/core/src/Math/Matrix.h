#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <stdexcept>

/**
 * @brief Alias for a sparse matrix using a map of coordinate pairs.
 * 
 * @tparam T Type of the row and column indices.
 */
template<typename T>
using SparseMatrix = std::map<std::pair<T, T>, double>;

/**
 * @brief Extracts the diagonal elements from a sparse matrix.
 * 
 * @tparam T Index type.
 * @param matrix The sparse matrix.
 * @return Vector of diagonal indices.
 */
template<typename T>
std::vector<T> getDiagonalElements(const SparseMatrix<T>& matrix) {
    std::vector<T> diagonalElements;
    for (const auto& element : matrix) {
        if (element.first.first == element.first.second) {
            diagonalElements.push_back(element.first.first);
        }
    }
    return diagonalElements;
}

/**
 * @brief Creates an identity matrix from a list of indices.
 * 
 * @tparam T Index type.
 * @param indices List of row/column indices.
 * @return Sparse identity matrix.
 */
template<typename T>
SparseMatrix<T> createIdentityMatrix(const std::vector<T>& indices) {
    SparseMatrix<T> identity;
    for (const T& index : indices) {
        identity[{index, index}] = 1.0;
    }
    return identity;
}

/**
 * @brief Retrieves an element from the sparse matrix.
 * 
 * @tparam T Index type.
 * @param matrix The sparse matrix.
 * @param row Row index.
 * @param col Column index.
 * @return Value of the element (0.0 if not present).
 */
template<typename T>
double getElement(const SparseMatrix<T>& matrix, const T& row, const T& col) {
    auto it = matrix.find({row, col});
    if (it != matrix.end()) {
        return it->second;
    }
    return 0.0;
}

/**
 * @brief Sets or updates an element in the sparse matrix.
 * 
 * @tparam T Index type.
 * @param matrix Sparse matrix to modify.
 * @param row Row index.
 * @param col Column index.
 * @param value Value to assign.
 */
template<typename T>
void setElement(SparseMatrix<T>& matrix, const T& row, const T& col, double value) {
    if (value != 0.0) {
        matrix[{row, col}] = value;
    } else {
        matrix.erase({row, col});
    }
}

/**
 * @brief Inverts a square sparse matrix using Gauss-Jordan elimination.
 * 
 * @tparam T Index type.
 * @param matrix Input square sparse matrix.
 * @param indices List of row/column indices (must match matrix size).
 * @return Inverted matrix.
 * @throws std::runtime_error if the matrix is singular.
 */
template<typename T>
SparseMatrix<T> invertMatrix(const SparseMatrix<T>& matrix, const std::vector<T>& indices) {
    SparseMatrix<T> augmentedMatrix = matrix;
    SparseMatrix<T> inverse = createIdentityMatrix(indices);

    for (size_t i = 0; i < indices.size(); ++i) {
        T pivot = indices[i];

        double pivotValue = getElement(augmentedMatrix, pivot, pivot);
        if (pivotValue == 0.0) {
            throw std::runtime_error("Matrix is singular and cannot be inverted.");
        }

        for (size_t k = 0; k < indices.size(); ++k) {
            T col = indices[k];
            setElement(augmentedMatrix, pivot, col, getElement(augmentedMatrix, pivot, col) / pivotValue);
            setElement(inverse, pivot, col, getElement(inverse, pivot, col) / pivotValue);
        }

        for (size_t j = 0; j < indices.size(); ++j) {
            if (i == j) continue;
            T row = indices[j];
            double factor = getElement(augmentedMatrix, row, pivot);
            for (size_t k = 0; k < indices.size(); ++k) {
                T col = indices[k];
                setElement(augmentedMatrix, row, col, getElement(augmentedMatrix, row, col) - factor * getElement(augmentedMatrix, pivot, col));
                setElement(inverse, row, col, getElement(inverse, row, col) - factor * getElement(inverse, pivot, col));
            }
        }
    }

    return inverse;
}

/**
 * @brief Prints the matrix in dense form using a list of indices.
 * 
 * @tparam T Index type.
 * @param matrix Sparse matrix.
 * @param indices List of indices to determine order.
 */
template<typename T>
void printMatrix(const SparseMatrix<T>& matrix, const std::vector<T>& indices) {
    for (const auto& i : indices) {
        for (const auto& j : indices) {
            std::cout << getElement(matrix, i, j) << " ";
        }
        std::cout << std::endl;
    }
}

#endif // __MATRIX_H__
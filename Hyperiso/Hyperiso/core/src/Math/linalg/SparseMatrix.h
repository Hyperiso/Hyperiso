#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

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
 * @brief Removes rows and columns that are entirely zero (i.e., not present in the sparse matrix).
 * 
 * @tparam T Index type.
 * @param matrix The input sparse matrix.
 * @param indices The input list of row/column indices.
 * @return A pair containing:
 *         - The cleaned sparse matrix.
 *         - The updated list of indices.
 */
template<typename T>
std::pair<SparseMatrix<T>, std::vector<T>> removeEmptyRowsAndCols(const SparseMatrix<T>& matrix, const std::vector<T>& indices) {
    std::map<T, int> rowCount;
    std::map<T, int> colCount;

    // Count non-zero elements per row and column
    for (const auto& [coord, value] : matrix) {
        if (value != 0.0) { // <-- Ignore zero entries explicitly present
            T row = coord.first;
            T col = coord.second;
            rowCount[row]++;
            colCount[col]++;
        }
    }

    // Keep only indices where row and column are not empty
    std::vector<T> cleanedIndices;
    for (const T& idx : indices) {
        if (rowCount[idx] > 0 || colCount[idx] > 0) {
            cleanedIndices.push_back(idx);
        }
    }

    // Filter matrix to keep only relevant rows/cols
    SparseMatrix<T> cleanedMatrix;
    for (const auto& [coord, value] : matrix) {
        if (value != 0.0) {
            T row = coord.first;
            T col = coord.second;
            if (std::find(cleanedIndices.begin(), cleanedIndices.end(), row) != cleanedIndices.end() &&
                std::find(cleanedIndices.begin(), cleanedIndices.end(), col) != cleanedIndices.end()) {
                cleanedMatrix[coord] = value;
            }
        }
    }

    return {cleanedMatrix, cleanedIndices};
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

/**
 * @brief Prints the matrix in a formatted grid with row and column headers.
 *
 * @tparam T Index type (must support std::ostream << operator).
 * @param matrix Sparse matrix.
 * @param indices List of row/column indices.
 */
template<typename T>
void customPrintMatrix(const SparseMatrix<T>& matrix, const std::vector<T>& indices) {
    // Determine maximum width for any label or value
    size_t maxLabelWidth = 0;
    size_t maxValueWidth = 0;

    std::ostringstream oss;
    for (const T& index : indices) {
        oss.str("");
        oss.clear();
        oss << index;
        maxLabelWidth = std::max(maxLabelWidth, oss.str().length());
    }

    for (const auto& [coord, value] : matrix) {
        oss.str("");
        oss.clear();
        oss << value;
        maxValueWidth = std::max(maxValueWidth, oss.str().length());
    }

    size_t cellWidth = std::max(maxLabelWidth, maxValueWidth) + 2;

    // Print top-left empty cell
    std::cout << std::setw(cellWidth) << " ";

    // Print column headers
    for (const auto& col : indices) {
        std::cout << std::setw(cellWidth) << col;
    }
    std::cout << "\n";

    // Print separator line
    std::cout << std::setw(cellWidth) << " ";
    for (size_t i = 0; i < indices.size(); ++i) {
        std::cout << std::setw(cellWidth) << std::string(cellWidth - 1, '-');
    }
    std::cout << "\n";

    // Print rows
    for (const auto& row : indices) {
        std::cout << std::setw(cellWidth) << row;
        for (const auto& col : indices) {
            double val = getElement(matrix, row, col);
            std::cout << std::setw(cellWidth) << val;
        }
        std::cout << "\n";
    }
}

template<typename T>
class SparseMatrixWrapper {
public:
    SparseMatrix<T> matrix;

    SparseMatrixWrapper() = default;

    double getElement(T row, T col) const {
        return ::getElement(matrix, row, col);
    }

    void setElement(T row, T col, double value) {
        ::setElement(matrix, row, col, value);
    }

    std::vector<T> getDiagonalElements() const {
        return ::getDiagonalElements(matrix);
    }

    void print(const std::vector<T>& indices) const {
        ::printMatrix(matrix, indices);
    }

    SparseMatrixWrapper<T> invert(const std::vector<T>& indices) const {
        SparseMatrixWrapper<T> result;
        result.matrix = ::invertMatrix(matrix, indices);
        return result;
    }
};


#endif // __SPARSEMATRIX_H__
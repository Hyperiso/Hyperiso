#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <stdexcept>
#include "Math.h"

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

template<typename T>
SparseMatrix<T> createIdentityMatrix(const std::vector<T>& indices) {
    SparseMatrix<T> identity;
    for (const T& index : indices) {
        identity[{index, index}] = 1.0;
    }
    return identity;
}

template<typename T>
double getElement(const SparseMatrix<T>& matrix, const T& row, const T& col) {
    auto it = matrix.find({row, col});
    if (it != matrix.end()) {
        return it->second;
    }
    return 0.0;
}

template<typename T>
void setElement(SparseMatrix<T>& matrix, const T& row, const T& col, double value) {
    if (value != 0.0) {
        matrix[{row, col}] = value;
    } else {
        matrix.erase({row, col});
    }
}

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

template<typename T>
void printMatrix(const SparseMatrix<T>& matrix, const std::vector<T>& indices) {
    for (const auto& i : indices) {
        for (const auto& j : indices) {
            std::cout << getElement(matrix, i, j) << " ";
        }
        std::cout << std::endl;
    }
}

// int main() {
//     SparseMatrix matrix = {
//         {{"1", "1"}, 1},
//         {{"1", "2"}, 3},
//         {{"2", "1"}, 3},
//         {{"2", "2"}, 1},
//         {{"3", "2"}, 5},
//         {{"3", "3"}, 1},
//         {{"2", "3"}, 5},
//         {{"d", "d"}, 1}
//     };

//     std::vector<std::string> indices = getDiagonalElements(matrix);

//     try {
//         SparseMatrix inverse = invertMatrix(matrix, indices);

//         printMatrix(inverse, indices);
//     } catch (const std::runtime_error& e) {
//         std::cerr << e.what() << std::endl;
//     }

//     return 0;
// }

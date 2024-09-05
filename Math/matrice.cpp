#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <stdexcept>
#include "Math.h"


std::vector<std::string> getDiagonalElements(const SparseMatrix& matrix) {
    std::vector<std::string> diagonalElements;
    for (const auto& element : matrix) {
        if (element.first.first == element.first.second) {
            diagonalElements.push_back(element.first.first);
        }
    }
    return diagonalElements;
}

SparseMatrix createIdentityMatrix(const std::vector<std::string>& indices) {
    SparseMatrix identity;
    for (const auto& index : indices) {
        identity[{index, index}] = 1.0;
    }
    return identity;
}

double getElement(const SparseMatrix& matrix, const std::string& row, const std::string& col) {
    auto it = matrix.find({row, col});
    if (it != matrix.end()) {
        return it->second;
    }
    return 0.0;
}

void setElement(SparseMatrix& matrix, const std::string& row, const std::string& col, double value) {
    if (value != 0.0) {
        matrix[{row, col}] = value;
    } else {
        matrix.erase({row, col});
    }
}

SparseMatrix invertMatrix(const SparseMatrix& matrix, const std::vector<std::string>& indices) {
    SparseMatrix augmentedMatrix = matrix;
    SparseMatrix inverse = createIdentityMatrix(indices);

    for (size_t i = 0; i < indices.size(); ++i) {
        std::string pivot = indices[i];

        double pivotValue = getElement(augmentedMatrix, pivot, pivot);
        if (pivotValue == 0.0) {
            throw std::runtime_error("Matrix is singular and cannot be inverted.");
        }

        for (size_t k = 0; k < indices.size(); ++k) {
            std::string col = indices[k];
            setElement(augmentedMatrix, pivot, col, getElement(augmentedMatrix, pivot, col) / pivotValue);
            setElement(inverse, pivot, col, getElement(inverse, pivot, col) / pivotValue);
        }

        for (size_t j = 0; j < indices.size(); ++j) {
            if (i == j) continue;
            std::string row = indices[j];
            double factor = getElement(augmentedMatrix, row, pivot);
            for (size_t k = 0; k < indices.size(); ++k) {
                std::string col = indices[k];
                setElement(augmentedMatrix, row, col, getElement(augmentedMatrix, row, col) - factor * getElement(augmentedMatrix, pivot, col));
                setElement(inverse, row, col, getElement(inverse, row, col) - factor * getElement(inverse, pivot, col));
            }
        }
    }

    return inverse;
}

void printMatrix(const SparseMatrix& matrix, const std::vector<std::string>& indices) {
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

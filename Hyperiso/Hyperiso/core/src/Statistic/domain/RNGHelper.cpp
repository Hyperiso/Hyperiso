#include "RNGHelper.h"

Matrix readMatrixFromStdin() {
    int n;
    if (!(std::cin >> n) || n <= 0) {
        throw std::runtime_error("Failed to read matrix size n.");
    }
    Matrix A(static_cast<size_t>(n), std::vector<double>(static_cast<size_t>(n)));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (!(std::cin >> A[i][j])) {
                throw std::runtime_error("Matrix lecture has failed.");
            }
        }
    }
    return A;
}

void printVector(const Vector& v) {
    std::cout << std::fixed << std::setprecision(15);
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) std::cout << " ";
        std::cout << v[i];
    }
    std::cout << "\n";
}

void printUsage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " [distribution=gaussian] [optional seed] < matrix.txt\n"
        << "  - The input matrix is read from standard input with the format:\n"
        << "      n\n"
        << "      r11 r12 ... r1n\n"
        << "      ...\n"
        << "      rn1 rn2 ... rnn\n"
        << "  - Supported distributions: gaussian | normal\n"
        << "Example:\n"
        << "  " << prog << " gaussian 12345 < my_corr.txt\n";
}

void printRealMatrix(const RealMatrix &A) {
    std::cout << std::scientific << std::setprecision(5);
    for (size_t i = 0; i < A.rows(); ++i) {
        std::cout << "[";
        for (size_t j = 0; j < A.cols(); j++) {
            std::cout << A.at(i, j);
            if (j != A.cols() - 1)
                std::cout << ", ";
        }
        std::cout << "]\n";
    }
}

#include "RNGHelper.h"

Matrix readMatrixFromStdin() {
    int n;
    if (!(std::cin >> n) || n <= 0) {
        throw std::runtime_error("Impossible to read n (matrix size).");
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
    std::cout << std::fixed << std::setprecision(6);
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) std::cout << " ";
        std::cout << v[i];
    }
    std::cout << "\n";
}

void printUsage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " [distribution=gaussian] [seed (optionnel)] < matrice.txt\n"
        << "  - La matrice d'entree est lue sur stdin au format:\n"
        << "      n\\n\n"
        << "      r11 r12 ... r1n\\n\n"
        << "      ...\\n"
        << "      rn1 rn2 ... rnn\\n\n"
        << "  - Distribution supportee: gaussian|normal\n"
        << "Exemple:\n"
        << "  " << prog << " gaussian 12345 < my_corr.txt\n";
}
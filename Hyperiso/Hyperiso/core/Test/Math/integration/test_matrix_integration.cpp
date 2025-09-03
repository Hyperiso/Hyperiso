// test_matrix_integration.cpp
#include <cassert>
#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <sstream>

#include "Matrix.h"

template <typename T>
static bool dapprox(double a, double b, double eps = 1e-12) {
    return std::abs(a - b) <= eps * (1.0 + std::abs(a) + std::abs(b));
}

template <typename T>
static SparseMatrix<T> multiply_dense(const SparseMatrix<T>& A, const SparseMatrix<T>& B, const std::vector<T>& idx) {
    SparseMatrix<T> C;
    for (const auto& i : idx) {
        for (const auto& j : idx) {
            double sum = 0.0;
            for (const auto& k : idx) {
                sum += getElement(A, i, k) * getElement(B, k, j);
            }
            if (sum != 0.0) C[{i,j}] = sum;
        }
    }
    return C;
}

template <typename T>
static bool is_identity(const SparseMatrix<T>& M, const std::vector<T>& idx, double eps = 1e-10) {
    for (const auto& i : idx) {
        for (const auto& j : idx) {
            double v = getElement(M, i, j);
            if (i == j) { if (!dapprox<T>(v, 1.0, eps)) return false; }
            else        { if (!dapprox<T>(v, 0.0, eps)) return false; }
        }
    }
    return true;
}

int main() {
    std::cout << "[integration] début des tests...\n";

    // 1) Inversion 3x3 + vérification A*A^{-1} ~ I avec indices permutés
    {
        // A = [[3,1,2],[0,1,4],[2,0,1]]
        SparseMatrix<int> A;
        setElement(A, 0,0, 3.0); setElement(A, 0,1, 1.0); setElement(A, 0,2, 2.0);
        setElement(A, 1,0, 0.0); setElement(A, 1,1, 1.0); setElement(A, 1,2, 4.0);
        setElement(A, 2,0, 2.0); setElement(A, 2,1, 0.0); setElement(A, 2,2, 1.0);
        // Indices dans un ordre non trié
        std::vector<int> idx{2,0,1};

        auto Ainv = invertMatrix(A, idx);
        auto I = multiply_dense(A, Ainv, idx);
        assert(is_identity(I, idx));
    }

    // 2) Nettoyage -> inversion
    {
        std::vector<int> idx{0,1,2,3,4};
        SparseMatrix<int> M;
        // Bloc 0..2 est utile, 3 et 4 sont vides
        setElement(M, 0,0, 10.0);
        setElement(M, 0,1, 2.0);
        setElement(M, 1,0, 3.0);
        setElement(M, 1,1, 9.0);
        setElement(M, 2,2, 5.0);

        // Nettoyage
        auto [Mc, idxc] = removeEmptyRowsAndCols(M, idx);
        assert(idxc == std::vector<int>({0,1,2}));

        // Inversion & check identité
        auto McInv = invertMatrix(Mc, idxc);
        auto I = multiply_dense(Mc, McInv, idxc);
        assert(is_identity(I, idxc));
    }

    // 3) Wrapper end-to-end
    {
        SparseMatrixWrapper<int> W;
        // A = [[2,1,0],[1,3,1],[0,1,2]]
        setElement(W.matrix, 0,0, 2.0); setElement(W.matrix, 0,1, 1.0);
        setElement(W.matrix, 1,0, 1.0); setElement(W.matrix, 1,1, 3.0); setElement(W.matrix, 1,2, 1.0);
        setElement(W.matrix, 2,1, 1.0); setElement(W.matrix, 2,2, 2.0);

        std::vector<int> idx{0,1,2};
        auto Winv = W.invert(idx);
        auto I = multiply_dense(W.matrix, Winv.matrix, idx);
        assert(is_identity(I, idx));

        // Impression (juste vérifier que ça écrit quelque chose)
        std::stringstream buf;
        std::streambuf* old = std::cout.rdbuf(buf.rdbuf());
        W.print(idx);
        std::cout.rdbuf(old);
        assert(!buf.str().empty());
    }

    std::cout << "All integration tests passed.\n";
    return 0;
}

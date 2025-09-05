#include <cassert>
#include <iostream>
#include <limits>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <set>
#include <complex>
#include <string>

#include "Matrix.h"

template <typename T>
static bool dapprox(T a, T b, double eps = 1e-12) {
    return std::abs(a - b) <= eps * (1.0 + std::abs(a) + std::abs(b));
}

template <typename T>
static bool mat_equal_approx(const SparseMatrix<T>& A, const SparseMatrix<T>& B, double eps = 1e-12) {
    if (A.size() != B.size()) return false;
    for (const auto& [k, va] : A) {
        auto it = B.find(k);
        if (it == B.end()) return false;
        if (!dapprox(va, it->second, eps)) return false;
    }
    return true;
}

int main() {
    std::cout << "[unit] début des tests matrix...\n";

    {
        std::vector<int> idx{0, 1, 3};
        auto I = createIdentityMatrix(idx);

        for (int i : idx) {
            assert(dapprox(getElement(I, i, i), 1.0));
            for (int j : idx) {
                if (i != j) assert(dapprox(getElement(I, i, j), 0.0));
            }
        }

        setElement(I, 0, 1, 2.5);
        assert(dapprox(getElement(I, 0, 1), 2.5));
        setElement(I, 0, 1, 0.0);
        assert(dapprox(getElement(I, 0, 1), 0.0));
        assert(I.find({0,1}) == I.end());
    }

    {
        std::vector<int> idx{2, 0, 1};
        auto I = createIdentityMatrix(idx);
        auto diag = getDiagonalElements(I);
        std::sort(diag.begin(), diag.end());
        std::sort(idx.begin(), idx.end());
        assert(diag == idx);
    }

    {
        SparseMatrix<int> A;
        setElement(A, 0, 0, 4.0);
        setElement(A, 0, 1, 7.0);
        setElement(A, 1, 0, 2.0);
        setElement(A, 1, 1, 6.0);
        std::vector<int> idx{0,1};

        auto Ainv = invertMatrix(A, idx);

        SparseMatrix<int> ref;
        setElement(ref, 0, 0, 0.6);
        setElement(ref, 0, 1, -0.7);
        setElement(ref, 1, 0, -0.2);
        setElement(ref, 1, 1, 0.4);

        assert(mat_equal_approx(Ainv, ref));
    }

    {
        SparseMatrix<int> S;
        setElement(S, 0, 0, 1.0);
        setElement(S, 0, 1, 2.0);
        setElement(S, 1, 0, 2.0);
        setElement(S, 1, 1, 4.0);
        std::vector<int> idx{0,1};

        bool thrown = false;
        try {
            auto Sinv = invertMatrix(S, idx);
            (void)Sinv;
        } catch (const std::runtime_error&) {
            thrown = true;
        }
        assert(thrown);
    }

    {
        std::vector<int> idx{0,1,2,3};

        SparseMatrix<int> M;
        setElement(M, 1, 1, 5.0);
        setElement(M, 1, 2, 1.0);
        setElement(M, 2, 1, -2.0);
        M[{3,3}] = 0.0;

        auto [Mc, idxc] = removeEmptyRowsAndCols(M, idx);

        std::vector<int> expectedIdx{1,2};
        assert(idxc == expectedIdx);

        assert(dapprox(getElement(Mc, 1, 1), 5.0));
        assert(dapprox(getElement(Mc, 1, 2), 1.0));
        assert(dapprox(getElement(Mc, 2, 1), -2.0));
        assert(dapprox(getElement(Mc, 2, 2), 0.0));
        assert(Mc.find({3,3}) == Mc.end());
    }

    {
        std::vector<int> idx{0,1};
        SparseMatrix<int> M;
        setElement(M, 0, 0, 1.0);
        setElement(M, 0, 1, 2.0);
        setElement(M, 1, 1, 3.0);

        std::stringstream buf;
        std::streambuf* old = std::cout.rdbuf(buf.rdbuf());
        printMatrix(M, idx);
        std::cout.rdbuf(old);

        const std::string expected = "1 2 \n0 3 \n";
        assert(buf.str() == expected);
    }

    {
        std::vector<int> idx{10, 20};
        SparseMatrix<int> M;
        setElement(M, 10, 10, 1.0);
        setElement(M, 10, 20, 2.0);
        setElement(M, 20, 10, 3.0);
        setElement(M, 20, 20, 4.0);

        std::stringstream buf;
        std::streambuf* old = std::cout.rdbuf(buf.rdbuf());
        customPrintMatrix(M, idx);
        std::cout.rdbuf(old);

        const std::string out = buf.str();
        assert(out.find("10") != std::string::npos);
        assert(out.find("20") != std::string::npos);
        assert(out.find("1")  != std::string::npos);
        assert(out.find("2")  != std::string::npos);
        assert(out.find("3")  != std::string::npos);
        assert(out.find("4")  != std::string::npos);
    }

    {
        SparseMatrixWrapper<int> W;
        W.setElement(0,0, 2.0);
        W.setElement(0,1, 1.0);
        W.setElement(1,0, 1.0);
        W.setElement(1,1, 2.0);

        assert(dapprox(W.getElement(0,0), 2.0));
        assert(dapprox(W.getElement(0,1), 1.0));

        auto diag = W.getDiagonalElements();
        std::sort(diag.begin(), diag.end());
        std::vector<int> expected{0,1};
        assert(diag == expected);

        auto Winv = W.invert({0,1});
        assert(dapprox(Winv.getElement(0,0), 2.0/3.0));
        assert(dapprox(Winv.getElement(0,1), -1.0/3.0));
        assert(dapprox(Winv.getElement(1,0), -1.0/3.0));
        assert(dapprox(Winv.getElement(1,1), 2.0/3.0));

        std::stringstream buf;
        std::streambuf* old = std::cout.rdbuf(buf.rdbuf());
        W.print({0,1});
        std::cout.rdbuf(old);
        assert(!buf.str().empty());
    }

    {
        using T = std::string;
        using namespace std::string_literals;

        std::vector<T> idx{"a"s, "b"s};
        SparseMatrix<T> A = createIdentityMatrix(idx);

        setElement(A, "a"s, "a"s, 2.0);
        setElement(A, "b"s, "b"s, 3.0);

        auto Ainverse = invertMatrix(A, idx);
        assert(dapprox(getElement(Ainverse, "a"s, "a"s), 0.5));
        assert(dapprox(getElement(Ainverse, "b"s, "b"s), 1.0/3.0));
        assert(dapprox(getElement(Ainverse, "a"s, "b"s), 0.0));
        assert(dapprox(getElement(Ainverse, "b"s, "a"s), 0.0));
    }

    std::cout << "All unit tests passed.\n";
    return 0;
}

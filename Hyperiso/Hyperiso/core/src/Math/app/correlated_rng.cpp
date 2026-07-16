#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

static constexpr double EPS_SYM = 1e-10;
static constexpr double EPS_DIAG = 1e-8;

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

struct IMarginalDistribution {
    virtual ~IMarginalDistribution() = default;

    virtual Vector sample(std::size_t n) = 0;
};

struct IDecomposition {
    virtual ~IDecomposition() = default;

    virtual Matrix factorize(const Matrix& R) = 0;
};

class CorrelationMatrixValidator {
public:
    void validate(const Matrix& R) const {
        const size_t n = R.size();
        if (n == 0) throw std::invalid_argument("Invalid matrix.");
        for (const auto& row : R) {
            if (row.size() != n) throw std::invalid_argument("Non squared matrix.");
        }

        //->check sym
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                if (std::fabs(R[i][j] - R[j][i]) > EPS_SYM) {
                    throw std::invalid_argument("Non symmetric matrix.");
                }
            }
        }
        //->check diag
        for (size_t i = 0; i < n; ++i) {
            if (std::fabs(R[i][i] - 1.0) > EPS_DIAG) {
                throw std::invalid_argument("Diagonal elements of the matrix needs to be ones");
            }
        }
        // positive-dev check with Cholesky
    }
};

// Cholesky 
class CholeskyDecomposition final : public IDecomposition {
public:
    Matrix factorize(const Matrix& R) override {
        const size_t n = R.size();
        Matrix L(n, std::vector<double>(n, 0.0));

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j <= i; ++j) {
                double sum = R[i][j];
                for (size_t k = 0; k < j; ++k) {
                    sum -= L[i][k] * L[j][k];
                }
                if (i == j) {
                    if (sum <= 0.0) {
                        throw std::invalid_argument(
                            "Matrix is not positive-definite for Cholesky decomposition.");
                    }
                    L[i][j] = std::sqrt(sum);
                } else {
                    L[i][j] = sum / L[j][j];
                }
            }
        }
        return L;
    }
};

class GaussianMarginal final : public IMarginalDistribution {
public:
    explicit GaussianMarginal(unsigned int seed = std::random_device{}())
        : eng_(seed), dist_(0.0, 1.0) {}

    Vector sample(std::size_t n) override {
        Vector z(n);
        for (std::size_t i = 0; i < n; ++i) z[i] = dist_(eng_);
        return z;
    }

private:
    std::mt19937 eng_;
    std::normal_distribution<double> dist_;
};

class DistributionFactory {
public:
    static std::unique_ptr<IMarginalDistribution> create(const std::string& name,
                                                 unsigned int seed = std::random_device{}()) {
        std::string lower = name;
        std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        });

        if (lower == "gaussian" || lower == "normal" || lower == "gauss") {
            return std::make_unique<GaussianMarginal>(seed);
        }

        throw std::invalid_argument("Unkwown distribution: " + name +
                                    " (try: gaussian|normal)");
    }
};

class JointDistribution {
public:
    JointDistribution(std::unique_ptr<IMarginalDistribution> dist,
                          std::unique_ptr<IDecomposition> decomp)
        : dist_(std::move(dist)), decomp_(std::move(decomp)) {}

    // y = L * z, z ~ i.i.d. (E=0, Var=1). Cov(y) = L L^T = R.
    Vector generate(const Matrix& correlation) const {
        CorrelationMatrixValidator().validate(correlation);
        Matrix L = decomp_->factorize(correlation);
        const std::size_t n = L.size();

        Vector z = dist_->sample(n);
        Vector y(n, 0.0);

        for (std::size_t i = 0; i < n; ++i) {
            double acc = 0.0;
            for (std::size_t k = 0; k <= i; ++k) {
                acc += L[i][k] * z[k];
            }
            y[i] = acc;
        }
        return y;
    }

private:
    std::unique_ptr<IMarginalDistribution> dist_;
    std::unique_ptr<IDecomposition> decomp_;
};

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

int main(int argc, char** argv) {
    try {
        std::string distName = "gaussian";
        unsigned int seed = std::random_device{}();

        if (argc >= 2) {
            std::string arg1 = argv[1];
            if (arg1 == "-h" || arg1 == "--help") {
                printUsage(argv[0]);
                return 0;
            }
            distName = arg1;
        }
        if (argc >= 3) {
            try {
                seed = static_cast<unsigned int>(std::stoul(argv[2]));
            } catch (...) {
                std::cerr << "Avertissement: seed invalide, utilisation d'un seed aleatoire.\n";
                seed = std::random_device{}();
            }
        }

        Matrix R = readMatrixFromStdin();

        auto dist = DistributionFactory::create(distName, seed);
        auto decomp = std::make_unique<CholeskyDecomposition>();

        JointDistribution generator(std::move(dist), std::move(decomp));
        Vector y = generator.generate(R);

        printVector(y);
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Erreur: " << ex.what() << "\n";
        printUsage(argv[0]);
        return 1;
    }
}

#include "JointDistribution.h"
#include "MarginalFactory.h"
#include "CopulaFactory.h"
#include "Math.h"

using Matrix = std::vector<std::vector<double>>;
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

        auto dist = MarginalFactory::create(MarginalType::GAUSSIAN, GaussianMarginalCfg(0., 1.), seed);
        // auto decomp = std::make_unique<CholeskyDecomposition>();
        auto copul = CopulaFactory::create(CopulaType::GAUSSIAN, GaussianCopulaConfig());

        std::vector<std::unique_ptr<IMarginalDistribution>> truc{};

        truc.emplace_back(std::move(dist));

        JointDistribution generator(std::move(truc), std::move(copul));
        Vector y = generator.sample();

        printVector(y);
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Erreur: " << ex.what() << "\n";
        printUsage(argv[0]);
        return 1;
    }
}

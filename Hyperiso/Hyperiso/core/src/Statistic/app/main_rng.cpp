#include "RandomVectorGenerator.h"
#include "DistributionFactory.h"

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

        auto dist = DistributionFactory::create(DistributionType::GAUSSIAN, seed);
        auto decomp = std::make_unique<CholeskyDecomposition>();

        RandomVectorGenerator generator(std::move(dist), std::move(decomp));
        Vector y = generator.generate(R);

        printVector(y);
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Erreur: " << ex.what() << "\n";
        printUsage(argv[0]);
        return 1;
    }
}

#include "JointDistribution.h"
#include "MarginalFactory.h"
#include "CopulaFactory.h"

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

        auto dist = DistributionFactory::create(MarginalType::GAUSSIAN, GaussianMarginalCfg(0., 1.), seed);
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

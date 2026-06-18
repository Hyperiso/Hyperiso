#include <iostream>
#include <string>

#include "Logger.h"
#include "WilsonHandler.h"
#include "ObservableHandler.h"
#include "StatisticHandler.h"

namespace {

void print_usage() {
    std::cout
        << "Hyperiso terminal interface\n\n"
        << "Usage:\n"
        << "  hyperiso-ui <module> <command> [options]\n\n"
        << "Modules:\n"
        << "  wilson      Builtin Wilson-coefficient summaries\n"
        << "  observable  Builtin observable summaries\n"
        << "  statistic   Statistic/dependency/uncertainty summaries\n\n"
        << "Common options:\n"
        << "  --model SM|THDM|MSSM|MARTY     Model, default SM\n"
        << "  --lha <path>                    LHA/FLHA input, default lha/si_input.flha\n"
        << "  --order LO|NLO|NNLO             QCD order, default NNLO\n"
        << "  --help                          Show help\n\n"
        << "Examples:\n"
        << "  hyperiso-ui wilson summary --groups BCoefficients --coeffs C7,C9,C10\n"
        << "  hyperiso-ui observable summary --observables BR_Bs__mu_mu,BR_B__Xs_gamma\n"
        << "  hyperiso-ui statistic summary --observables BR_Bs__mu_mu --draws 200 --progress\n";
}

} // namespace

int main(int argc, char* argv[]) {
    Logger::getInstance()->setEnabled(true);

    if (argc < 2 || std::string(argv[1]) == "--help" || std::string(argv[1]) == "help") {
        print_usage();
        return argc < 2 ? 1 : 0;
    }

    const std::string module = argv[1];
    try {
        if (module == "wilson") return handleWilsonOptions(argc - 1, argv + 1);
        if (module == "observable") return handleObservableOptions(argc - 1, argv + 1);
        if (module == "statistic") return handleStatisticOptions(argc - 1, argv + 1);

        std::cerr << "Unknown module: " << module << "\n\n";
        print_usage();
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}

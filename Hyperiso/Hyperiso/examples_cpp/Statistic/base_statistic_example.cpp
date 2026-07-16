#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "HyperisoMaster.h"
#include "Include.h"
#include "Logger.h"
#include "ObservableInterface.h"
#include "StatisticInterface.h"

namespace {

void print_fit_result(const FitResultWithMaps& fit) {
    std::cout << "fit_ok = " << fit.fit_ok << "\n";
    std::cout << "ell_hat = " << fit.ell_hat << "\n";

    std::cout << "p_hat:\n";
    for (const auto& [pid, value] : fit.p_hat) {
        std::cout << "  " << pid << " = " << value << "\n";
    }

    std::cout << "p_hat_std:\n";
    for (const auto& [pid, value] : fit.p_hat_std) {
        std::cout << "  " << pid << " = " << value << "\n";
    }
}

} // namespace

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    HyperisoConfig config;
    config.model = Model::SM;

    HyperisoMaster hyp;
    hyp.init("lha/si_input.flha", config);

    auto oi = std::make_shared<ObservableInterface>();

    // Add all observables from a decay. The last argument adds dependencies, useful for statistics workflows.
    oi->add_observables(Decays::B__l_l, QCDOrder::NNLO, true);

    // Configuration for the statistical interface.
    // A lot of options can be fine-tuned but we advise the user only to change few of them.
    StatisticConfig sc;

    // The most useful parameter to change is the number of draws for the Monte Carlo generator.
    sc.MC_draws = 100;

    // Main interface for uncertainty calculation, MLE and making 2D contours.
    StatisticInterface si(sc, oi);

    // This API calculates symmetric and asymmetric uncertainties of the observables given in the ObservableInterface.
    // It calculates the skew as well to know if the distribution is symmetric or not.
    const auto result_u = si.compute_uncertainties();

    std::cout << "\nUncertainties:\n";
    for (const auto& [obs, summary] : result_u) {
        std::cout << obs.str() << " -> " << summary << "\n";
    }

    // One can also fit over one or multiple parameters through the following API.
    // It returns the best values for the p_specs and the nuisances.
    const std::vector<ParamId> p_specs = {
        ParamId(ParameterType::FLAVOR, "FCONST", LhaID(511, 1)),
        ParamId(ParameterType::FLAVOR, "FCONST", LhaID(531, 1))
    };

    const auto result_mle = si.compute_MLE(p_specs);

    std::cout << "\nMLE:\n";
    print_fit_result(result_mle);

    // ContourOptions can be tuned if you want another contour algorithm or a different resolution.
    ContourOptions co;

    const auto result_contour = si.compute_confidence_contour(
        ParamId(ParameterType::FLAVOR, "FCONST", LhaID(511, 1)),
        ParamId(ParameterType::FLAVOR, "FCONST", LhaID(531, 1)),
        1.0,
        {-0.5, 0.5, -0.5, 0.5},
        co
    );

    std::cout << "\nContour success: " << result_contour.success
              << ", level: " << result_contour.level
              << ", paths: " << result_contour.paths.size() << "\n";

    return 0;
}

#include <iostream>
#include <memory>
#include <vector>

#include "HyperisoMaster.h"
#include "Include.h"
#include "Logger.h"
#include "ObservableInterface.h"
#include "StatisticInterface.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    // Dev example: prepare a likelihood once, then evaluate a small 2D scan.
    HyperisoConfig config;
    config.model = Model::SM;

    HyperisoMaster hyp;
    hyp.init("lha/si_input.flha", config);

    auto oi = std::make_shared<ObservableInterface>();
    oi->add_observable(Observables::BR_BS_MUMU, QCDOrder::NNLO, true)
       .add_observable(Observables::BR_BD_MUMU, QCDOrder::NNLO, true);

    StatisticConfig sc;
    sc.MC_draws = 100;
    sc.MLE_max_iter = 500;

    StatisticInterface si(sc, oi);

    const ParamId f_bd(ParameterType::FLAVOR, "FCONST", LhaID(511, 1));
    const ParamId f_bs(ParameterType::FLAVOR, "FCONST", LhaID(531, 1));
    const std::vector<ParamId> p_specs = {f_bd, f_bs};

    // Build and store the likelihood objects needed for scan calls.
    si.prepare_likelihood_for_scan(p_specs);

    // Evaluate a rectangular grid around the current point.
    const auto grid = si.scan_likelihood_around_current_point(
        f_bd,
        f_bs,
        0.03,   // half width for f_Bd
        0.03,   // half width for f_Bs
        7,      // grid points on x
        7       // grid points on y
    );

    std::cout << "Scan for " << grid.x_param << " and " << grid.y_param << "\n";
    std::cout << "center = (" << grid.x_center << ", " << grid.y_center << ")\n";
    std::cout << "grid size = " << grid.nx << " x " << grid.ny << "\n";

    std::cout << "\nFirst points:\n";
    std::size_t printed = 0;
    for (const auto& p : grid.points) {
        if (printed++ >= 10) break;
        std::cout << "  x=" << p.x
                  << ", y=" << p.y
                  << ", nll=" << p.nll
                  << ", delta_nll=" << p.delta_nll << "\n";
    }

    // You can save the grid as CSV for plotting in Python, ROOT or gnuplot.
    si.save_likelihood_scan_csv("likelihood_scan.csv", grid);
    std::cout << "Wrote likelihood_scan.csv\n";

    return 0;
}

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "StatisticManager.h"
#include "ObservableInterface.h"
#include "ObservableInterfaceProxy.h"
#include "StatCorrelationProxy.h"
#include "StatParameterProxy.h"
#include "StatParamSourcesProxy.h"
#include "StatDependencyPruner.h"
#include "FitAbstraction.h"
#include "NuisanceReader.h"
#include "DefaultNuisancePathsProvider.h"

namespace {

using Path = std::vector<std::pair<double, double>>;

void print_fit_result(const FitResultWithMaps& fit) {
    std::cout << std::setprecision(17);
    std::cout << "fit_ok  = " << fit.fit_ok << "\n";
    std::cout << "ell_hat = " << fit.ell_hat << "\n";

    std::cout << "\np_hat:\n";
    for (const auto& [pid, val] : fit.p_hat) {
        std::cout << "  " << fit_app::to_string_any(pid) << " = " << val << "\n";
    }

    std::cout << "\neta_hat:\n";
    for (const auto& [pid, val] : fit.eta_hat) {
        std::cout << "  " << fit_app::to_string_any(pid) << " = " << val << "\n";
    }

    std::cout << "\np_hat_std:\n";
    for (const auto& [pid, val] : fit.p_hat_std) {
        std::cout << "  " << fit_app::to_string_any(pid) << " = " << val << "\n";
    }

    std::cout << "\np_correlations:\n";
    for (const auto& [pi, row] : fit.p_correlations) {
        for (const auto& [pj, corr] : row) {
            std::cout << "  corr(" << fit_app::to_string_any(pi)
                      << ", " << fit_app::to_string_any(pj)
                      << ") = " << corr << "\n";
        }
    }
}

void save_bestfit_csv(const std::string& path, const FitResultWithMaps& fit) {
    std::ofstream out(path);
    out << "name,value,error\n";

    for (const auto& [pid, val] : fit.p_hat) {
        double err = 0.0;
        auto it = fit.p_hat_std.find(pid);
        if (it != fit.p_hat_std.end()) {
            err = it->second;
        }

        out << fit_app::to_string_any(pid) << ","
            << std::setprecision(17) << val << ","
            << err << "\n";
    }
}

void save_contour_csv(const std::string& path,
                      const std::string& xname,
                      const std::string& yname,
                      const std::set<Path>& contour68,
                      const std::set<Path>& contour95) {
    std::ofstream out(path);
    out << "# x=" << xname << "\n";
    out << "# y=" << yname << "\n";
    out << "cl,path_id,x,y\n";

    std::size_t path_id = 0;
    for (const auto& path_pts : contour68) {
        for (const auto& p : path_pts) {
            out << "0.683," << path_id << ","
                << std::setprecision(17) << p.first << ","
                << p.second << "\n";
        }
        ++path_id;
    }

    path_id = 0;
    for (const auto& path_pts : contour95) {
        for (const auto& p : path_pts) {
            out << "0.95," << path_id << ","
                << std::setprecision(17) << p.first << ","
                << p.second << "\n";
        }
        ++path_id;
    }
}

} 

int main() {
    using namespace fit_app;
    // Logger::getInstance()->setLevel(Logger::LogLevel::VERBOSE);
    HyperisoMaster hyp;
    HyperisoConfig config_hyp;
    config_hyp.model = Model::SM;
    hyp.init("lha/si_input.flha", config_hyp);

    auto oint = std::make_shared<ObservableInterface>();
    oint->add_observable(ObservableMapper::to_id(Observables::BR_BS_MUMU_UNTAG), QCDOrder::LO, true)
        .add_observable(ObservableMapper::to_id(Observables::BR_BD_MUMU), QCDOrder::LO, true);

    std::shared_ptr<IStatParamOptimizerProxy> spop = std::make_shared<StatParamOptimizerProxy>();

    auto model = std::make_shared<ObservableInterfaceProxy>(oint, spop);

    StatisticConfig config;
    config.MC_draws = 100;
    config.MLE_max_iter = 120000;
    config.MLE_tol = 0.2;

    std::vector<ParamId> p_specs = {
        ParamId{ParameterType::FLAVOR, "FCONST", {511, 1}},
        ParamId{ParameterType::FLAVOR, "FCONST", {531, 1}}
    };

    std::shared_ptr<INuisancePathsProvider> npp = std::make_shared<DefaultNuisancePathsProvider>();

    StatisticManager stat(
        config,
        model,
        std::make_shared<StatCorrelationProxy>(),
        std::make_shared<StatParameterProxy>(),
        std::make_shared<StatParamSourcesProxy>(),
        std::make_shared<StatDependencyPruner>(),
        std::make_shared<NuisanceReader>(npp),
        spop
    );

    auto t0 = std::chrono::steady_clock::now();
    // auto unc = stat.compute_uncertainties();
    auto t1 = std::chrono::steady_clock::now();
    std::cout << "Uncertainty pass done in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
              << " ms\n";
    // std::cout << "Number of summarized observables = " << unc.size() << "\n";

    auto t2 = std::chrono::steady_clock::now();
    FitResultWithMaps fit = stat.compute_MLE(p_specs);
    auto t3 = std::chrono::steady_clock::now();

    std::cout << "\nMLE done in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count()
              << " ms\n\n";

    if (!fit.fit_ok) {
        std::cerr << "[ERROR] MLE fit failed.\n";
        return 5;
    }

    print_fit_result(fit);
    save_bestfit_csv("bestfit.csv", fit);
    std::cout << "[INFO] Wrote bestfit.csv\n";

    if (p_specs.size() == 2) {
        const ParamId p1 = p_specs[0];
        const ParamId p2 = p_specs[1];

        std::array<double, 4> bounds = {0.05, 0.35, 0.05, 0.35};

        //   68.3%  -> delta_chi2 = 2.30 -> delta_nll = 1.15 -> z = sqrt(2.30)
        //   95%    -> delta_chi2 = 5.99 -> delta_nll = 2.995 -> z = sqrt(5.99)
        const double z68_2d = std::sqrt(2.30);
        const double z95_2d = std::sqrt(5.99);

        ContourOptions  opt;
        
        opt.fallback_contour_method = ContourAlgorithm::AMS;
        opt.primary_contour_method = ContourAlgorithm::MINUIT;
        opt.profile_backend= ProfileBackend::LAPLACE_NUISANCE;
        
        auto c68 = stat.confidence_contour(p1, p2, z68_2d, bounds, opt);

        std::cout << "first contour done" << std::endl; 

        auto t4 = std::chrono::steady_clock::now();
        auto c95 = stat.confidence_contour(p1, p2, z95_2d, bounds, opt);
        auto t5 = std::chrono::steady_clock::now();

        std::cout << "\nContour done in "
                << std::chrono::duration_cast<std::chrono::milliseconds>(t5 - t4).count()
                << " ms\n\n";

        

        std::cout << "[INFO] contour 68% paths = " << c68.level << "\n";
        std::cout << "[INFO] contour 95% paths = " << c95.level << "\n";

        save_contour_csv(
            "contours.csv",
            fit_app::to_string_any(p1),
            fit_app::to_string_any(p2),
            c68.paths,
            c95.paths
        );
        std::cout << "[INFO] Wrote contours.csv\n";
    }

    return 0;
}
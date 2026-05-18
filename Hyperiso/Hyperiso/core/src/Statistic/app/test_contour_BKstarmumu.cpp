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
#include "ObservableInterfaceAdapter2.h"
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
    // oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_1_CPV_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_1_CPV_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_1_CPV_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
    //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {6, 8.68}}, QCDOrder::NNLO, true)
    //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_1_CPV_B0__KSTAR0_MU_MU), {6, 8.68}}, QCDOrder::NNLO, true)
    //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), {6, 8.68}}, QCDOrder::NNLO, true)
    //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), {6, 8.68}}, QCDOrder::NNLO, true)
    //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), {6, 8.68}}, QCDOrder::NNLO, true)
    //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), {6, 8.68}}, QCDOrder::NNLO, true)
    //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), {6, 8.68}}, QCDOrder::NNLO, true)
    //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), {6, 8.68}}, QCDOrder::NNLO, true)
    //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {10.09, 12.86}}, QCDOrder::NNLO, true)
    //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_1_CPV_B0__KSTAR0_MU_MU), {10.09, 12.86}}, QCDOrder::NNLO, true)
    //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), {10.09, 12.86}}, QCDOrder::NNLO, true)
    //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), {10.09, 12.86}}, QCDOrder::NNLO, true)
    //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), {10.09, 12.86}}, QCDOrder::NNLO, true)
    //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), {10.09, 12.86}}, QCDOrder::NNLO, true)
    //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), {10.09, 12.86}}, QCDOrder::NNLO, true)
    //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), {10.09, 12.86}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_1_CPV_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
    //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true);

        const std::vector<Observables> ang_obs = {
            // Observables::DGAMMA_DQ2_B0__KSTAR0_MU_MU,
            Observables::F_L_B0__KSTAR0_MU_MU,
            Observables::P_1_B0__KSTAR0_MU_MU,
            Observables::P_2_B0__KSTAR0_MU_MU,
            Observables::P_3_B0__KSTAR0_MU_MU,
            Observables::P_PRIME_4_B0__KSTAR0_MU_MU,
            Observables::P_PRIME_5_B0__KSTAR0_MU_MU,
            Observables::P_PRIME_6_B0__KSTAR0_MU_MU,
            Observables::P_PRIME_8_B0__KSTAR0_MU_MU,
            Observables::S_1C_B0__KSTAR0_MU_MU,
            Observables::S_2S_B0__KSTAR0_MU_MU
        };

        auto add_ang_bin = [&](double q2min, double q2max) {
            for (auto obs : ang_obs) {
                oint->add_observable(
                    BinnedObservableId{ObservableMapper::to_id(obs), {q2min, q2max}},
                    QCDOrder::NNLO,
                    true
                );
            }
        };

        // CMS, Table 1, en excluant le bin problématique [6, 8.68]
        add_ang_bin(1.1, 2.0);
        add_ang_bin(2.0, 4.3);
        add_ang_bin(4.3, 6.0);
        // add_ang_bin(6.0, 8.68); // à exclure pour reproduire Table 1
        // add_ang_bin(10.09, 12.86);
        add_ang_bin(14.18, 16.0);
        // add_ang_bin(0.06, 0.98);
        // add_ang_bin(1.1, 2.5);
        // add_ang_bin(2.5, 4.0);
        // add_ang_bin(4.0, 6.0);
        // // add_ang_bin(6.0, 8.0); // à exclure
        // add_ang_bin(15.0, 17.0);
        // add_ang_bin(17.0, 19.0);

        // // // // LHCb2025 config 2, en excluant [6, 8]
        // add_ang_bin(0.06, 0.98);
        // add_ang_bin(1.1, 2.5);
        // add_ang_bin(2.5, 4.0);
        // add_ang_bin(4.0, 6.0);
        // // // add_ang_bin(6.0, 8.0); // à exclure
        // add_ang_bin(15.0, 17.0);
        // add_ang_bin(17.0, 19.0);
        // oint->add_observable(
        //             BinnedObservableId{ObservableMapper::to_id(Observables::S_2S_B0__KSTAR0_MU_MU), {0.06, 0.98}},
        //             QCDOrder::NNLO,
        //             true
        //         );
        // oint->add_observable(
        //             BinnedObservableId{ObservableMapper::to_id(Observables::S_6C_B0__KSTAR0_MU_MU), {0.06, 0.98}},
        //             QCDOrder::NNLO,
        //             true
        //         );


    BKstarllConfig cfg;
    cfg.ff_src = BV_FF_Src::GRvDV;
    oint->set_decay_config(Decays::B__Kstar_l_l, cfg);
    oint->set_bkstarll_threads(20);
    std::shared_ptr<IStatParamOptimizerProxy> spop = std::make_shared<StatParamOptimizerProxy>();
    auto model = std::make_shared<ObservableInterfaceAdapterObs>(oint, spop);
    

    StatisticConfig config;
    config.MLE_max_iter = 120000;
    config.MLE_tol = 0.1;
    config.MLE_trace_first_evals  = true;
    config.MLE_trace_max_evals  = 20;


    const std::string had_bsm_block =
        GroupMapper::str(WGroup::B, ScaleType::HADRONIC, WilsonBasis::B_STANDARD)
        + "__BSM_INTERMEDIATE";

    std::vector<ParamId> p_specs = {
        ParamId{ParameterType::WILSON, had_bsm_block, WCoefMapper::flha_full(WCoef::C9, QCDOrder::LO, ContributionType::BSM)},
        ParamId{ParameterType::WILSON, had_bsm_block, WCoefMapper::flha_full(WCoef::C10, QCDOrder::LO, ContributionType::BSM)}
    };

    // std::vector<ParamId> p_specs = {
    //     // ParamId{ParameterType::WILSON, had_bsm_block, WCoefMapper::flha_full(WCoef::C9, QCDOrder::LO, ContributionType::BSM)},
    //     ParamId{ParameterType::WILSON, had_bsm_block, WCoefMapper::flha_full(WCoef::C10, QCDOrder::LO, ContributionType::BSM)}
    // };

    LOG_INFO("Creating StatisticManager.");

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
    // std::set<std::string> exp = {"CMS", "LHCb2020", "LHCb2025c2"};
    // std::set<std::string> exp = {"LHCb2020"};
    // std::set<std::string> exp = {"LHCb2025c2"};
    std::set<std::string> exp = {"CMS"};
    // std::set<std::string> exp = {"CMS", "LHCb2020"};
    stat.select_experiments(exp);
    
    LOG_INFO("StatisticManager created.");

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

        std::array<double, 4> bounds = {-8, 8, -8, 8};

        auto trace = std::make_shared<std::ofstream>("contour_trace.csv");
        (*trace) << "type,level,path_id,point_id,x,y,n_paths,n_points,elapsed_s,message\n";

        ContourOptions  opt;
        opt.primary_contour_method = ContourAlgorithm::AMS;
        opt.fallback_contour_method = ContourAlgorithm::MINUIT;
        opt.profile_backend = ProfileBackend::LAPLACE_NUISANCE;
        opt.resolution = 20;

        opt.on_progress = [trace](const ContourProgressEvent& ev) {
        if (!trace || !(*trace)) return;

        const char* type_str = "";
        switch (ev.type) {
            case ContourProgressEventType::Started:      type_str = "started"; break;
            case ContourProgressEventType::PathPoint:    type_str = "point"; break;
            case ContourProgressEventType::PathFinished: type_str = "path_finished"; break;
            case ContourProgressEventType::Finished:     type_str = "finished"; break;
            case ContourProgressEventType::Failed:       type_str = "failed"; break;
        }

        (*trace)
            << type_str << ","
            << ev.level << ","
            << ev.path_id << ","
            << ev.point_id << ","
            << ev.x << ","
            << ev.y << ","
            << ev.n_paths << ","
            << ev.n_points << ","
            << ev.elapsed_seconds << ","
            << "\"" << ev.message << "\"\n";

        trace->flush();
    };
        auto c68 = stat.confidence_contour(p1, p2, 1, bounds, opt);
        std::cout << "finish first contour" << std::endl;

        auto t4 = std::chrono::steady_clock::now();
        auto c95 = stat.confidence_contour(p1, p2, 2, bounds, opt);
        auto t5 = std::chrono::steady_clock::now();

        std::cout << "\nContour done in "
                << std::chrono::duration_cast<std::chrono::milliseconds>(t5 - t4).count()
                << " ms\n\n";

        std::cout << "[INFO] contour 68% paths = " << c68.level << "\n";
        std::cout << "[INFO] contour 95% paths = " << c95.level << "\n";

        save_contour_csv(
            "contours_BKsmumu_ang.csv",
            fit_app::to_string_any(p1),
            fit_app::to_string_any(p2),
            c68.paths,
            c95.paths
        );
        std::cout << "[INFO] Wrote contours.csv\n";
    }

    return 0;
}
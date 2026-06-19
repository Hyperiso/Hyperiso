// #include <chrono>
// #include <fstream>
// #include <iomanip>
// #include <iostream>
// #include <map>
// #include <memory>
// #include <set>
// #include <string>
// #include <utility>
// #include <vector>

// #include "StatisticManager.h"
// #include "ObservableInterface.h"
// #include "ObservableInterfaceAdapter2.h"
// #include "StatCorrelationProxy.h"
// #include "StatParameterProxy.h"
// #include "StatParamSourcesProxy.h"
// #include "StatDependencyPruner.h"
// #include "FitAbstraction.h"
// #include "NuisanceReader.h"
// #include "DefaultNuisancePathsProvider.h"

// namespace {

// using Path = std::vector<std::pair<double, double>>;

// void print_fit_result(const FitResultWithMaps& fit) {
//     std::cout << std::setprecision(17);
//     std::cout << "fit_ok  = " << fit.fit_ok << "\n";
//     std::cout << "ell_hat = " << fit.ell_hat << "\n";

//     std::cout << "\np_hat:\n";
//     for (const auto& [pid, val] : fit.p_hat) {
//         std::cout << "  " << fit_app::to_string_any(pid) << " = " << val << "\n";
//     }

//     std::cout << "\neta_hat:\n";
//     for (const auto& [pid, val] : fit.eta_hat) {
//         std::cout << "  " << fit_app::to_string_any(pid) << " = " << val << "\n";
//     }

//     std::cout << "\np_hat_std:\n";
//     for (const auto& [pid, val] : fit.p_hat_std) {
//         std::cout << "  " << fit_app::to_string_any(pid) << " = " << val << "\n";
//     }

//     std::cout << "\np_correlations:\n";
//     for (const auto& [pi, row] : fit.p_correlations) {
//         for (const auto& [pj, corr] : row) {
//             std::cout << "  corr(" << fit_app::to_string_any(pi)
//                       << ", " << fit_app::to_string_any(pj)
//                       << ") = " << corr << "\n";
//         }
//     }
// }

// void save_bestfit_csv(const std::string& path, const FitResultWithMaps& fit) {
//     std::ofstream out(path);
//     out << "name,value,error\n";

//     for (const auto& [pid, val] : fit.p_hat) {
//         double err = 0.0;
//         auto it = fit.p_hat_std.find(pid);
//         if (it != fit.p_hat_std.end()) {
//             err = it->second;
//         }

//         out << fit_app::to_string_any(pid) << ","
//             << std::setprecision(17) << val << ","
//             << err << "\n";
//     }
// }

// void save_contour_csv(const std::string& path,
//                       const std::string& xname,
//                       const std::string& yname,
//                       const std::set<Path>& contour68,
//                       const std::set<Path>& contour95) {
//     std::ofstream out(path);
//     out << "# x=" << xname << "\n";
//     out << "# y=" << yname << "\n";
//     out << "cl,path_id,x,y\n";

//     std::size_t path_id = 0;
//     for (const auto& path_pts : contour68) {
//         for (const auto& p : path_pts) {
//             out << "0.683," << path_id << ","
//                 << std::setprecision(17) << p.first << ","
//                 << p.second << "\n";
//         }
//         ++path_id;
//     }

//     path_id = 0;
//     for (const auto& path_pts : contour95) {
//         for (const auto& p : path_pts) {
//             out << "0.95," << path_id << ","
//                 << std::setprecision(17) << p.first << ","
//                 << p.second << "\n";
//         }
//         ++path_id;
//     }
// }

// } 

// int main() {
//     using namespace fit_app;
//     // Logger::getInstance()->setLevel(Logger::LogLevel::VERBOSE);
//     HyperisoMaster hyp;
//     HyperisoConfig config_hyp;
//     config_hyp.model = Model::SM;
//     hyp.init("lha/si_input.flha", config_hyp);

//     auto oint = std::make_shared<ObservableInterface>();
//     // oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_1_CPV_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_1_CPV_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_1_CPV_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
//     //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {6, 8.68}}, QCDOrder::NNLO, true)
//     //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_1_CPV_B0__KSTAR0_MU_MU), {6, 8.68}}, QCDOrder::NNLO, true)
//     //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), {6, 8.68}}, QCDOrder::NNLO, true)
//     //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), {6, 8.68}}, QCDOrder::NNLO, true)
//     //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), {6, 8.68}}, QCDOrder::NNLO, true)
//     //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), {6, 8.68}}, QCDOrder::NNLO, true)
//     //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), {6, 8.68}}, QCDOrder::NNLO, true)
//     //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), {6, 8.68}}, QCDOrder::NNLO, true)
//     //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {10.09, 12.86}}, QCDOrder::NNLO, true)
//     //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_1_CPV_B0__KSTAR0_MU_MU), {10.09, 12.86}}, QCDOrder::NNLO, true)
//     //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), {10.09, 12.86}}, QCDOrder::NNLO, true)
//     //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), {10.09, 12.86}}, QCDOrder::NNLO, true)
//     //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), {10.09, 12.86}}, QCDOrder::NNLO, true)
//     //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), {10.09, 12.86}}, QCDOrder::NNLO, true)
//     //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), {10.09, 12.86}}, QCDOrder::NNLO, true)
//     //     // .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), {10.09, 12.86}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_1_CPV_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
//     //     .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true);

//         // const std::vector<Observables> ang_obs = {
//         //     // Observables::DGAMMA_DQ2_B0__KSTAR0_MU_MU,
//         //     Observables::F_L_B0__KSTAR0_MU_MU,
//         //     Observables::P_1_B0__KSTAR0_MU_MU,
//         //     Observables::P_2_B0__KSTAR0_MU_MU,
//         //     Observables::P_3_B0__KSTAR0_MU_MU,
//         //     Observables::P_PRIME_4_B0__KSTAR0_MU_MU,
//         //     Observables::P_PRIME_5_B0__KSTAR0_MU_MU,
//         //     Observables::P_PRIME_6_B0__KSTAR0_MU_MU,
//         //     Observables::P_PRIME_8_B0__KSTAR0_MU_MU,
//         //     Observables::S_1C_B0__KSTAR0_MU_MU,
//         //     // Observables::S_2S_B0__KSTAR0_MU_MU,
//         //     // Observables::S_6C_B0__KSTAR0_MU_MU
//         // };
//         BKstarllConfig cfg;
//         cfg.ff_src = BV_FF_Src::GRvDV;
//         oint->set_decay_config(Decays::B__Kstar_l_l, cfg);
//         oint->set_bkstarll_threads(16);
//         const std::vector<Observables> ang_obs_lhcb2025c2 = {
//             Observables::DGAMMA_DQ2_B0__KSTAR0_MU_MU,
//             Observables::F_L_B0__KSTAR0_MU_MU,
//             Observables::S_1C_B0__KSTAR0_MU_MU,
//             Observables::P_1_B0__KSTAR0_MU_MU,
//             Observables::P_2_B0__KSTAR0_MU_MU,
//             Observables::P_3_B0__KSTAR0_MU_MU,
//             Observables::P_PRIME_4_B0__KSTAR0_MU_MU,
//             Observables::P_PRIME_5_B0__KSTAR0_MU_MU,
//             Observables::P_PRIME_6_B0__KSTAR0_MU_MU,
//             Observables::P_PRIME_8_B0__KSTAR0_MU_MU
//         };

//         const std::vector<Observables> ang_obs_cms = {
//             // Observables::DGAMMA_DQ2_B0__KSTAR0_MU_MU,
//             Observables::F_L_B0__KSTAR0_MU_MU,
//             Observables::P_1_B0__KSTAR0_MU_MU,
//             Observables::P_2_B0__KSTAR0_MU_MU,
//             Observables::P_3_B0__KSTAR0_MU_MU,
//             Observables::P_PRIME_4_B0__KSTAR0_MU_MU,
//             Observables::P_PRIME_5_B0__KSTAR0_MU_MU,
//             Observables::P_PRIME_6_B0__KSTAR0_MU_MU,
//             Observables::P_PRIME_8_B0__KSTAR0_MU_MU
//         };

//         // ang_obs = DecayMapper::get_observables(Decays::B__Kstar_l_l);
//         // auto add_ang_bin = [&](double q2min, double q2max) {
//         //     for (auto obs : ang_obs) {
//         //         oint->add_observable(
//         //             BinnedObservableId{ObservableMapper::to_id(obs), {q2min, q2max}},
//         //             QCDOrder::NNLO,
//         //             true
//         //         );
//         //     }
//         // };
        

//         using O = Observables;

//         // Évite les doublons exacts obs/bin, utile quand DEFAULT, Belle ou CMS ont parfois le même bin.
//         std::set<Observables> seen_unbinned;
//         std::set<std::tuple<Observables, double, double>> seen_binned;

//         auto add_unbinned = [&](Observables obs) {
//             if (!seen_unbinned.insert(obs).second) return;

//             oint->add_observable(
//                 ObservableMapper::to_id(obs),
//                 QCDOrder::NNLO,
//                 true
//             );
//         };

//         auto add_bin = [&](Observables obs, double q2min, double q2max) {
//             auto key = std::make_tuple(obs, q2min, q2max);
//             if (!seen_binned.insert(key).second) return;

//             oint->add_observable(
//                 BinnedObservableId{ObservableMapper::to_id(obs), {q2min, q2max}},
//                 QCDOrder::NNLO,
//                 true
//             );
//         };

//         auto add_bins = [&](std::initializer_list<Observables> obs_list,
//                             std::initializer_list<std::pair<double, double>> bins) {
//             for (auto obs : obs_list) {
//                 for (auto [q2min, q2max] : bins) {
//                     add_bin(obs, q2min, q2max);
//                 }
//             }
//         };

//         // ============================================================
//         // Unbinned / radiatifs / leptoniques
//         // ============================================================

//         add_unbinned(O::IA_B__KSTAR_GAMMA);       // AI_BKstargamma dans ton tableau
//         add_unbinned(O::BR_B_XS_GAMMA);
//         add_unbinned(O::BR_BS_MUMU_UNTAG);
//         add_unbinned(O::BR_BS_EE_UNTAG);
//         add_unbinned(O::BR_B0__KSTAR0_GAMMA);
//         add_unbinned(O::BR_B__KSTAR_GAMMA);

//         // ============================================================
//         // B -> Xs ll
//         // ============================================================

//         // add_bins(
//         //     {
//         //         O::BR_B__Xs_mu_mu,
//         //         O::BR_B__Xs_e_e
//         //     },
//         //     {
//         //         {1.0, 6.0},
//         //         {14.2, 22.0}
//         //     }
//         // );

//         add_unbinned(O::BR_B__Xs_mu_mu);
//         add_unbinned(O::BR_B__Xs_e_e);

//         // ============================================================
//         // B -> K* mu mu
//         // ============================================================

//         add_bins(
//             {O::DBR_DQ2_B__KSTAR_MU_MU},
//             {
//                 {1.1, 6.0},
//                 {15.0, 19.0}
//             }
//         );

//         add_bins(
//             {
//                 O::F_L_B__KSTAR_MU_MU,
//                 O::A_FB_B__KSTAR_MU_MU,
//                 O::S_3_B__KSTAR_MU_MU,
//                 O::S_4_B__KSTAR_MU_MU,
//                 O::S_5_B__KSTAR_MU_MU,
//                 O::S_7_B__KSTAR_MU_MU,
//                 O::S_8_B__KSTAR_MU_MU,
//                 O::S_9_B__KSTAR_MU_MU
//             },
//             {
//                 {0.1, 0.98},
//                 {1.1, 2.5},
//                 {2.5, 4.0},
//                 {4.0, 6.0},
//                 {6.0, 8.0},
//                 {15.0, 17.0},
//                 {17.0, 19.0}
//             }
//         );

//         // ============================================================
//         // Ratios R - 1
//         // ============================================================

//         add_bins(
//             {O::R_1_B0__KSTAR0_L_L},
//             {
//                 {0.1, 1.1},
//                 {1.1, 6.0},
//                 {0.045, 1.1},  // Belle
//                 {15.0, 19.0}   // Belle
//             }
//         );

//         add_bins(
//             {O::R_1_B__K_L_L},
//             {
//                 {0.1, 1.1},
//                 {1.1, 6.0},
//                 {1.0, 6.0}     // Belle
//             }
//         );

//         add_bin(O::R_1_B__KSTAR_L_L, 0.045, 6.0);

//         // MANQUANT DANS GeneralEnum.h :
//         // R-1_B0K0ll_1.1_6 devrait correspondre à un enum du style
//         // O::R_1_B0__K0_L_L, mais cet enum n'existe pas dans le fichier uploadé.
//         // add_bin(O::R_1_B0__K0_L_L, 1.1, 6.0);

//         // ============================================================
//         // B -> K mu mu
//         // ============================================================

//         add_bins(
//             {O::DGAMMA_DQ2_B0__K0_MU_MU},
//             {
//                 {1.1, 6.0},
//                 {15.0, 22.0}
//             }
//         );

//         add_bins(
//             {O::DGAMMA_DQ2_B__K_MU_MU},
//             {
//                 {1.1, 6.0},
//                 {15.0, 22.0}
//             }
//         );

//         add_bins(
//             {O::F_H_B__K_MU_MU},
//             {
//                 {1.1, 6.0},
//                 {15.0, 22.0}
//             }
//         );

//         // ============================================================
//         // Bs -> phi mu mu
//         // ============================================================

//         add_bins(
//             {O::DGAMMA_DQ2_BS__PHI_MU_MU},
//             {
//                 {0.1, 0.98},
//                 {1.1, 2.5},
//                 {2.5, 4.0},
//                 {4.0, 6.0},
//                 {6.0, 8.0},
//                 {15.0, 19.0}
//             }
//         );

//         add_bins(
//             {
//                 O::F_L_BS_PHI_MU_MU,
//                 O::S_3_BS_PHI_MU_MU,
//                 O::S_4_BS_PHI_MU_MU,
//                 O::S_7_BS_PHI_MU_MU
//             },
//             {
//                 {0.1, 0.98},
//                 {1.1, 4.0},
//                 {4.0, 6.0},
//                 {6.0, 8.0},
//                 {15.0, 18.9}
//             }
//         );

//         // ============================================================
//         // Lambda_b -> Lambda mu mu
//         // ============================================================

//         add_bins(
//             {
//                 O::DGAMMA_DQ2_LAMBDA_B__LAMBDA_MU_MU,
//                 O::A_FB_L_LAMBDA_B__LAMBDA_MU_MU,
//                 O::A_FB_H_LAMBDA_B__LAMBDA_MU_MU,
//                 O::A_FB_LH_LAMBDA_B__LAMBDA_MU_MU,
//                 O::F_L_LAMBDA_B__LAMBDA_MU_MU
//             },
//             {
//                 {15.0, 20.0}
//             }
//         );

//         // ============================================================
//         // B0 -> K*0 ee
//         // ============================================================

//         add_bin(O::DGAMMA_DQ2_B0__KSTAR0_E_E, 0.0009, 1.0);

//         add_bins(
//             {
//                 O::F_L_B0__KSTAR0_E_E,
//                 O::A_T_RE_B0__KSTAR0_E_E,
//                 O::A_T_2_B0__KSTAR0_E_E
//             },
//             {
//                 {0.0008, 0.257}
//             }
//         );

//         // Belle
//         add_bin(O::A_T_2_B0__KSTAR0_E_E, 0.0008, 1.12);

//         add_bins(
//             {
//                 O::F_L_B0__KSTAR0_E_E,
//                 O::P_1_B0__KSTAR0_E_E,
//                 O::P_2_B0__KSTAR0_E_E,
//                 O::P_3_B0__KSTAR0_E_E,
//                 O::P_PRIME_4_B0__KSTAR0_E_E,
//                 O::P_PRIME_5_B0__KSTAR0_E_E,
//                 O::P_PRIME_6_B0__KSTAR0_E_E,
//                 O::P_PRIME_8_B0__KSTAR0_E_E
//             },
//             {
//                 {1.1, 6.0}
//             }
//         );

//         // ============================================================
//         // Bs -> phi ee
//         // ============================================================

//         add_bins(
//             {O::DGAMMA_DQ2_BS__PHI_E_E},
//             {
//                 {0.1, 1.1},
//                 {1.1, 6.0},
//                 {15.0, 19.0}
//             }
//         );

//         add_bins(
//             {O::R_1_BS__PHI_L_L},
//             {
//                 {0.1, 1.1},
//                 {1.1, 6.0},
//                 {15.0, 19.0}
//             }
//         );

//         add_bins(
//             {
//                 O::F_L_BS_PHI_E_E,
//                 O::A_T_2_BS_PHI_E_E
//             },
//             {
//                 {0.0009, 0.2615}
//             }
//         );

//         // ============================================================
//         // CMS : B -> K mu mu
//         // ============================================================

//         add_bin(O::F_H_B__K_MU_MU, 1.0, 6.0);

//         add_bins(
//             {O::DGAMMA_DQ2_B__K_MU_MU},
//             {
//                 {0.1, 0.98},
//                 {1.1, 2.0},
//                 {2.0, 3.0},
//                 {3.0, 4.0},
//                 {4.0, 5.0},
//                 {5.0, 6.0},
//                 {6.0, 7.0},
//                 {7.0, 8.0},
//                 {14.82, 16.0},
//                 {16.0, 17.0},
//                 {17.0, 18.0},
//                 {18.0, 19.24},
//                 {19.24, 22.9}
//             }
//         );

//         add_bin(O::R_1_B__K_L_L, 1.1, 6.0);

//         // ============================================================
//         // CMS : B0 -> K*0 mu mu
//         // ============================================================

//         add_bins(
//             {
//                 O::F_L_B0__KSTAR0_MU_MU,
//                 O::P_1_B0__KSTAR0_MU_MU,
//                 O::P_2_B0__KSTAR0_MU_MU,
//                 O::P_3_B0__KSTAR0_MU_MU,
//                 O::P_PRIME_4_B0__KSTAR0_MU_MU,
//                 O::P_PRIME_5_B0__KSTAR0_MU_MU,
//                 O::P_PRIME_6_B0__KSTAR0_MU_MU,
//                 O::P_PRIME_8_B0__KSTAR0_MU_MU
//             },
//             {
//                 {1.1, 2.0},
//                 {2.0, 4.3},
//                 {4.3, 6.0},
//                 {6.0, 8.68},
//                 {14.18, 16.0}
//             }
//         );

//         // ============================================================
//         // B0 -> K*0 mu mu : bin bas [0.06, 0.98]
//         // ============================================================

//         add_bins(
//             {
//                 O::F_L_B0__KSTAR0_MU_MU,
//                 O::S_2S_B0__KSTAR0_MU_MU,
//                 O::S_1C_B0__KSTAR0_MU_MU,
//                 O::P_1_B0__KSTAR0_MU_MU,
//                 O::P_2_B0__KSTAR0_MU_MU,
//                 O::P_3_B0__KSTAR0_MU_MU,
//                 O::P_PRIME_4_B0__KSTAR0_MU_MU,
//                 O::P_PRIME_5_B0__KSTAR0_MU_MU,
//                 O::P_PRIME_6_B0__KSTAR0_MU_MU,
//                 O::P_PRIME_8_B0__KSTAR0_MU_MU,
//                 O::S_6C_B0__KSTAR0_MU_MU,
//                 O::DGAMMA_DQ2_B0__KSTAR0_MU_MU
//             },
//             {
//                 {0.06, 0.98}
//             }
//         );

//         // ============================================================
//         // LHCb2025c2 : B0 -> K*0 mu mu
//         // ============================================================

//         add_bins(
//             {
//                 O::F_L_B0__KSTAR0_MU_MU,
//                 O::S_1C_B0__KSTAR0_MU_MU,
//                 O::P_1_B0__KSTAR0_MU_MU,
//                 O::P_2_B0__KSTAR0_MU_MU,
//                 O::P_3_B0__KSTAR0_MU_MU,
//                 O::P_PRIME_4_B0__KSTAR0_MU_MU,
//                 O::P_PRIME_5_B0__KSTAR0_MU_MU,
//                 O::P_PRIME_6_B0__KSTAR0_MU_MU,
//                 O::P_PRIME_8_B0__KSTAR0_MU_MU,
//                 O::DGAMMA_DQ2_B0__KSTAR0_MU_MU
//             },
//             {
//                 {1.1, 2.5},
//                 {2.5, 4.0},
//                 {4.0, 6.0},
//                 {6.0, 8.0},
//                 {15.0, 17.0},
//                 {17.0, 19.0}
//             }
//         );

//         // add_ang_bin(ang_obs_lhcb2025c2, 0.06, 0.98);
//         // add_ang_bin(ang_obs_lhcb2025c2, 1.1, 2.5);
//         // add_ang_bin(ang_obs_lhcb2025c2, 2.5, 4.0);
//         // add_ang_bin(ang_obs_lhcb2025c2, 4.0, 6.0);
//         // add_ang_bin(ang_obs_lhcb2025c2, 15.0, 17.0);
//         // add_ang_bin(ang_obs_lhcb2025c2, 17.0, 19.0);

//         // oint->add_observable(
//         //     BinnedObservableId{
//         //         ObservableMapper::to_id(Observables::S_2S_B0__KSTAR0_MU_MU),
//         //         {0.06, 0.98}
//         //     },
//         //     QCDOrder::NNLO,
//         //     true
//         // );

//         // oint->add_observable(
//         //     BinnedObservableId{
//         //         ObservableMapper::to_id(Observables::S_6C_B0__KSTAR0_MU_MU),
//         //         {0.06, 0.98}
//         //     },
//         //     QCDOrder::NNLO,
//         //     true
//         // );

//         // add_ang_bin(ang_obs_cms, 1.1, 2.0);
//         // add_ang_bin(ang_obs_cms, 2.0, 4.3);
//         // add_ang_bin(ang_obs_cms, 4.3, 6.0);
//         // add_ang_bin(ang_obs_cms, 14.18, 16.0);

//         // // // CMS, Table 1, en excluant le bin problématique [6, 8.68]
//         // // add_ang_bin(1.1, 2.0);
//         // // add_ang_bin(2.0, 4.3);
//         // // add_ang_bin(4.3, 6.0);
//         // // // add_ang_bin(6.0, 8.68); // à exclure pour reproduire Table 1
//         // // // add_ang_bin(10.09, 12.86);
//         // // add_ang_bin(14.18, 16.0);
//         // add_ang_bin(0.06, 0.98);
//         // // add_ang_bin(0.1, 0.98);
//         // // add_ang_bin(1.1, 1.6);
//         // add_ang_bin(1.1, 2.5);
//         // add_ang_bin(2.5, 4.0);
//         // add_ang_bin(4.0, 6.0);
//         // // add_ang_bin(6.0, 8.0); // à exclure
//         // // add_ang_bin(11.0, 12.5);
//         // add_ang_bin(15.0, 17.0);
//         // add_ang_bin(17.0, 19.0);
//         // // add_ang_bin(15.0, 19.0);

//         // // // // LHCb2025 config 2, en excluant [6, 8]
//         // add_ang_bin(0.06, 0.98);
//         // add_ang_bin(1.1, 2.5);
//         // add_ang_bin(2.5, 4.0);
//         // add_ang_bin(4.0, 6.0);
//         // // // add_ang_bin(6.0, 8.0); // à exclure
//         // add_ang_bin(15.0, 17.0);
//         // add_ang_bin(17.0, 19.0);
//         // oint->add_observable(
//         //             BinnedObservableId{ObservableMapper::to_id(Observables::S_2S_B0__KSTAR0_MU_MU), {0.06, 0.98}},
//         //             QCDOrder::NNLO,
//         //             true
//         //         );
//         // oint->add_observable(
//         //             BinnedObservableId{ObservableMapper::to_id(Observables::S_6C_B0__KSTAR0_MU_MU), {0.06, 0.98}},
//         //             QCDOrder::NNLO,
//         //             true
//         //         );


    
//     std::shared_ptr<IStatParamOptimizerProxy> spop = std::make_shared<StatParamOptimizerProxy>();
//     auto model = std::make_shared<ObservableInterfaceAdapterObs>(oint, spop);
    

//     StatisticConfig config;
//     config.advanced.MLE_max_iter = 120000;
//     config.advanced.MLE_tol = 0.05;
//     config.advanced.MLE_trace_first_evals  = true;
//     config.advanced.MLE_trace_max_evals  = 20;
//     config.advanced.likelihood_mode = StatisticLikelihoodMode::CHI2_MC_COVARIANCE;
//     config.MC_draws = 1000;
//     const std::string had_bsm_block =
//         GroupMapper::str(WGroup::B, ScaleType::HADRONIC, WilsonBasis::B_STANDARD)
//         + "__BSM_INTERMEDIATE";

//     std::vector<ParamId> p_specs = {
//         ParamId{ParameterType::WILSON, had_bsm_block, WCoefMapper::flha_full(WCoef::C9, QCDOrder::LO, ContributionType::BSM)},
//         ParamId{ParameterType::WILSON, had_bsm_block, WCoefMapper::flha_full(WCoef::C10, QCDOrder::LO, ContributionType::BSM)}
//     };

//     // std::vector<ParamId> p_specs = {
//     //     // ParamId{ParameterType::WILSON, had_bsm_block, WCoefMapper::flha_full(WCoef::C9, QCDOrder::LO, ContributionType::BSM)},
//     //     ParamId{ParameterType::WILSON, had_bsm_block, WCoefMapper::flha_full(WCoef::C10, QCDOrder::LO, ContributionType::BSM)}
//     // };

//     LOG_INFO("Creating StatisticManager.");

//     std::shared_ptr<INuisancePathsProvider> npp = std::make_shared<DefaultNuisancePathsProvider>();

//     StatisticManager stat(
//         config,
//         model,
//         std::make_shared<StatCorrelationProxy>(),
//         std::make_shared<StatParameterProxy>(),
//         std::make_shared<StatParamSourcesProxy>(),
//         std::make_shared<StatDependencyPruner>(),
//         std::make_shared<NuisanceReader>(npp),
//         spop
//     );
//     // std::set<std::string> exp = {"CMS", "LHCb2025c2"};
//     std::set<std::string> exp = {"DEFAULT", "Belle", "CMS", "LHCb2025c2"};
//     // std::set<std::string> exp = {"LHCb2020"};
//     // std::set<std::string> exp = {"LHCb2025c2"};
//     // std::set<std::string> exp = {"CMS"};
//     // std::set<std::string> exp = {"CMS", "LHCb2020"};
//     stat.select_experiments(exp);
    
//     LOG_INFO("StatisticManager created.");

//     auto t2 = std::chrono::steady_clock::now();
//     FitResultWithMaps fit = stat.compute_MLE(p_specs);
//     auto t3 = std::chrono::steady_clock::now();

//     std::cout << "\nMLE done in "
//               << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count()
//               << " ms\n\n";

//     if (!fit.fit_ok) {
//         std::cerr << "[ERROR] MLE fit failed.\n";
//         return 5;
//     }

//     print_fit_result(fit);
//     save_bestfit_csv("bestfit.csv", fit);
//     std::cout << "[INFO] Wrote bestfit.csv\n";

//     if (p_specs.size() == 2) {
//         const ParamId p1 = p_specs[0];
//         const ParamId p2 = p_specs[1];

//         // std::array<double, 4> bounds = {-20, 20, -20, 20};
//         std::array<double, 4> bounds = {
//             -3.0, 2.0,
//             -2.5, 2.5
//         };

//         auto trace = std::make_shared<std::ofstream>("contour_trace.csv");
//         (*trace) << "type,level,path_id,point_id,x,y,n_paths,n_points,elapsed_s,message\n";

//         ContourOptions  opt;
//         opt.primary_contour_method = ContourAlgorithm::MINUIT;
//         opt.fallback_contour_method = ContourAlgorithm::AMS;
//         opt.profile_backend = ProfileBackend::MINUIT;
//         opt.resolution = 10;

//         opt.on_progress = [trace](const ContourProgressEvent& ev) {
//         if (!trace || !(*trace)) return;

//         const char* type_str = "";
//         switch (ev.type) {
//             case ContourProgressEventType::Started:      type_str = "started"; break;
//             case ContourProgressEventType::PathPoint:    type_str = "point"; break;
//             case ContourProgressEventType::PathFinished: type_str = "path_finished"; break;
//             case ContourProgressEventType::Finished:     type_str = "finished"; break;
//             case ContourProgressEventType::Failed:       type_str = "failed"; break;
//         }

//         (*trace)
//             << type_str << ","
//             << ev.level << ","
//             << ev.path_id << ","
//             << ev.point_id << ","
//             << ev.x << ","
//             << ev.y << ","
//             << ev.n_paths << ","
//             << ev.n_points << ","
//             << ev.elapsed_seconds << ","
//             << "\"" << ev.message << "\"\n";

//         trace->flush();
//     };
//         auto t4 = std::chrono::steady_clock::now();
//         auto c95 = stat.confidence_contour(p1, p2, 2, bounds, opt);
//         auto t5 = std::chrono::steady_clock::now();
        
//         auto c68 = stat.confidence_contour(p1, p2, 1, bounds, opt);
//         std::cout << "finish first contour" << std::endl;


//         std::cout << "\nContour done in "
//                 << std::chrono::duration_cast<std::chrono::milliseconds>(t5 - t4).count()
//                 << " ms\n\n";

//         std::cout << "[INFO] contour 68% paths = " << c68.level << "\n";
//         std::cout << "[INFO] contour 95% paths = " << c95.level << "\n";

//         save_contour_csv(
//             "contours_BKsmumu_ang.csv",
//             fit_app::to_string_any(p1),
//             fit_app::to_string_any(p2),
//             c68.paths,
//             c95.paths
//         );
//         std::cout << "[INFO] Wrote contours.csv\n";
//     }

//     return 0;
// }

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
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
#include "BlockProxy.h"

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

    

    BlockProxy().log_block(ParameterType::DECAY, "B_Xsll");

    auto oint = std::make_shared<ObservableInterface>();
        // ============================================================
        // Observables: STRICTEMENT ceux du fichier "Texte collé(7).txt"
        // ============================================================

        BKstarllConfig cfg_BKs;
        cfg_BKs.ff_src = BV_FF_Src::GRvDV;
        oint->set_decay_config(Decays::B__Kstar_l_l, cfg_BKs);
        oint->set_bkstarll_threads(24);

        BKstarGammaConfig cfg_BKsgamma;
        cfg_BKsgamma.ff_src = BV_FF_Src::GRvDV;
        oint->set_decay_config(Decays::B__Kstar_gamma, cfg_BKsgamma);

        BsPhiConfig cfg_BsPhi;
        cfg_BsPhi.ff_src = BV_FF_Src::GRvDV;
        oint->set_decay_config(Decays::Bs__phi_l_l, cfg_BsPhi);
        oint->set_bsphi_threads(24);

        BKllConfig cfg_BK;
        cfg_BK.ff_src = BP_FF_Src::GKvD_SR_LAT;
        oint->set_decay_config(Decays::B__K_l_l, cfg_BK);
        oint->set_bkll_threads(24);

        using O = Observables;
        constexpr bool kAddDeps = false;

        std::set<Observables> seen_unbinned;
        std::set<std::tuple<Observables, double, double>> seen_binned;

        auto add_unbinned = [&](Observables obs) {
            if (!seen_unbinned.insert(obs).second) return;
            oint->add_observable(
                ObservableMapper::to_id(obs),
                QCDOrder::NNLO,
                kAddDeps
            );
        };

        auto add_bin = [&](Observables obs, double q2min, double q2max) {
            auto key = std::make_tuple(obs, q2min, q2max);
            if (!seen_binned.insert(key).second) return;
            oint->add_observable(
                BinnedObservableId{ObservableMapper::to_id(obs), {q2min, q2max}},
                QCDOrder::NNLO,
                kAddDeps
            );
        };

        add_unbinned(O::IA_B__KSTAR_GAMMA); // 001 AI_BKstargamma
        add_unbinned(O::BR_B_XS_GAMMA); // 002 BR_BXsgamma
        add_unbinned(O::BR_BS_MUMU_UNTAG); // 003 BRuntag_Bsmumu
        add_unbinned(O::BR_BS_EE_UNTAG); // 004 BRuntag_Bsee
        add_bin(O::BR_B__Xs_mu_mu, 1, 6); // 005 BR_BXsmumu_1_6
        add_bin(O::BR_B__Xs_mu_mu, 14.2, 22); // 006 BR_BXsmumu_14.2_22
        add_bin(O::BR_B__Xs_e_e, 1, 6); // 007 BR_BXsee_1_6
        add_bin(O::BR_B__Xs_e_e, 14.2, 22); // 008 BR_BXsee_14.2_22
        add_unbinned(O::BR_B0__KSTAR0_GAMMA); // 009 BR_B0Kstar0gamma
        add_unbinned(O::BR_B__KSTAR_GAMMA); // 010 BR_BKstargamma
        add_bin(O::DGAMMA_DQ2_B__KSTAR_MU_MU, 1.1, 6); // 011 dGamma/dq2_BKstarmumu_1.1_6
        add_bin(O::DBR_DQ2_B__KSTAR_MU_MU, 15, 19); // 012 dGamma/dq2_BKstarmumu_15_19
        add_bin(O::R_1_B0__KSTAR0_L_L, 0.1, 1.1); // 013 R-1_B0Kstar0ll_0.1_1.1
        add_bin(O::R_1_B0__KSTAR0_L_L, 1.1, 6); // 014 R-1_B0Kstar0ll_1.1_6
        add_bin(O::R_1_B0__KSTAR0_L_L, 0.045, 1.1); // 015 R-1_B0Kstar0ll_0.045_1.1_Belle
        add_bin(O::R_1_B0__KSTAR0_L_L, 1.1, 6); // 016 R-1_B0Kstar0ll_1.1_6_Belle
        add_bin(O::R_1_B0__KSTAR0_L_L, 15, 19); // 017 R-1_B0Kstar0ll_15_19_Belle
        add_bin(O::DBR_DQ2_B0__K0_MU_MU, 1.1, 6); // 018 dGamma/dq2_B0K0mumu_1.1_6
        add_bin(O::DBR_DQ2_B0__K0_MU_MU, 15, 22); // 019 dGamma/dq2_B0K0mumu_15_22
        add_bin(O::DBR_DQ2_B__K_MU_MU, 1.1, 6); // 020 dGamma/dq2_BKmumu_1.1_6
        add_bin(O::F_H_B__K_MU_MU, 1.1, 6); // 021 FH_BKmumu_1.1_6
        add_bin(O::DBR_DQ2_B__K_MU_MU, 15, 22); // 022 dGamma/dq2_BKmumu_15_22
        add_bin(O::F_H_B__K_MU_MU, 15, 22); // 023 FH_BKmumu_15_22
        add_bin(O::R_1_B__K_L_L, 0.1, 1.1); // 024 R-1_BKll_0.1_1.1
        add_bin(O::R_1_B__K_L_L, 1.1, 6); // 025 R-1_BKll_1.1_6
        add_bin(O::DBR_DQ2_BS__PHI_MU_MU, 0.1, 0.98); // 026 dGamma/dq2_Bsphimumu_0.1_0.98
        add_bin(O::F_L_BS_PHI_MU_MU, 0.1, 0.98); // 027 FL_Bsphimumu_0.1_0.98
        add_bin(O::S_3_BS_PHI_MU_MU, 0.1, 0.98); // 028 S3_Bsphimumu_0.1_0.98
        add_bin(O::S_4_BS_PHI_MU_MU, 0.1, 0.98); // 029 S4_Bsphimumu_0.1_0.98
        add_bin(O::S_7_BS_PHI_MU_MU, 0.1, 0.98); // 030 S7_Bsphimumu_0.1_0.98
        add_bin(O::DBR_DQ2_BS__PHI_MU_MU, 1.1, 2.5); // 031 dGamma/dq2_Bsphimumu_1.1_2.5
        add_bin(O::DBR_DQ2_BS__PHI_MU_MU, 2.5, 4); // 032 dGamma/dq2_Bsphimumu_2.5_4
        add_bin(O::F_L_BS_PHI_MU_MU, 1.1, 4); // 033 FL_Bsphimumu_1.1_4
        add_bin(O::S_3_BS_PHI_MU_MU, 1.1, 4); // 034 S3_Bsphimumu_1.1_4
        add_bin(O::S_4_BS_PHI_MU_MU, 1.1, 4); // 035 S4_Bsphimumu_1.1_4
        add_bin(O::S_7_BS_PHI_MU_MU, 1.1, 4); // 036 S7_Bsphimumu_1.1_4
        add_bin(O::DBR_DQ2_BS__PHI_MU_MU, 4, 6); // 037 dGamma/dq2_Bsphimumu_4_6
        add_bin(O::F_L_BS_PHI_MU_MU, 4, 6); // 038 FL_Bsphimumu_4_6
        add_bin(O::S_3_BS_PHI_MU_MU, 4, 6); // 039 S3_Bsphimumu_4_6
        add_bin(O::S_4_BS_PHI_MU_MU, 4, 6); // 040 S4_Bsphimumu_4_6
        add_bin(O::S_7_BS_PHI_MU_MU, 4, 6); // 041 S7_Bsphimumu_4_6
        add_bin(O::DBR_DQ2_BS__PHI_MU_MU, 15, 19); // 042 dGamma/dq2_Bsphimumu_15_19 // [2026-05-24_01] [WARN] Rejected MC nuisance sample 1 while trying to fill accepted sample 1 of 1000 : MC prediction contains non-finite observable
        add_bin(O::F_L_BS_PHI_MU_MU, 15, 18.9); // 043 FL_Bsphimumu_15_18.9
        add_bin(O::S_3_BS_PHI_MU_MU, 15, 18.9); // 044 S3_Bsphimumu_15_18.9
        add_bin(O::S_4_BS_PHI_MU_MU, 15, 18.9); // 045 S4_Bsphimumu_15_18.9
        add_bin(O::S_7_BS_PHI_MU_MU, 15, 18.9); // 046 S7_Bsphimumu_15_18.9
        add_bin(O::DBR_DQ2_LAMBDA_B__LAMBDA_MU_MU, 15, 20); // 047 dGamma/dq2_LambdabLambdamumu_15_20
        add_bin(O::A_FB_L_LAMBDA_B__LAMBDA_MU_MU, 15, 20); // 048 AlFB_LambdabLambdamumu_15_20
        add_bin(O::A_FB_H_LAMBDA_B__LAMBDA_MU_MU, 15, 20); // 049 AhFB_LambdabLambdamumu_15_20
        add_bin(O::A_FB_LH_LAMBDA_B__LAMBDA_MU_MU, 15, 20); // 050 AlhFB_LambdabLambdamumu_15_20
        add_bin(O::F_L_LAMBDA_B__LAMBDA_MU_MU, 15, 20); // 051 FL_LambdabLambdamumu_15_20
        add_bin(O::F_L_B__KSTAR_MU_MU, 0.1, 0.98); // 052 FL_BKstarmumu_0.1_0.98
        add_bin(O::A_FB_B__KSTAR_MU_MU, 0.1, 0.98); // 053 AFB_BKstarmumu_0.1_0.98
        add_bin(O::S_3_B__KSTAR_MU_MU, 0.1, 0.98); // 054 S3_BKstarmumu_0.1_0.98
        add_bin(O::S_4_B__KSTAR_MU_MU, 0.1, 0.98); // 055 S4_BKstarmumu_0.1_0.98
        add_bin(O::S_5_B__KSTAR_MU_MU, 0.1, 0.98); // 056 S5_BKstarmumu_0.1_0.98
        add_bin(O::S_7_B__KSTAR_MU_MU, 0.1, 0.98); // 057 S7_BKstarmumu_0.1_0.98
        add_bin(O::S_8_B__KSTAR_MU_MU, 0.1, 0.98); // 058 S8_BKstarmumu_0.1_0.98
        add_bin(O::S_9_B__KSTAR_MU_MU, 0.1, 0.98); // 059 S9_BKstarmumu_0.1_0.98
        add_bin(O::F_L_B__KSTAR_MU_MU, 1.1, 2.5); // 060 FL_BKstarmumu_1.1_2.5
        add_bin(O::A_FB_B__KSTAR_MU_MU, 1.1, 2.5); // 061 AFB_BKstarmumu_1.1_2.5
        add_bin(O::S_3_B__KSTAR_MU_MU, 1.1, 2.5); // 062 S3_BKstarmumu_1.1_2.5
        add_bin(O::S_4_B__KSTAR_MU_MU, 1.1, 2.5); // 063 S4_BKstarmumu_1.1_2.5
        add_bin(O::S_5_B__KSTAR_MU_MU, 1.1, 2.5); // 064 S5_BKstarmumu_1.1_2.5
        add_bin(O::S_7_B__KSTAR_MU_MU, 1.1, 2.5); // 065 S7_BKstarmumu_1.1_2.5
        add_bin(O::S_8_B__KSTAR_MU_MU, 1.1, 2.5); // 066 S8_BKstarmumu_1.1_2.5
        add_bin(O::S_9_B__KSTAR_MU_MU, 1.1, 2.5); // 067 S9_BKstarmumu_1.1_2.5
        add_bin(O::F_L_B__KSTAR_MU_MU, 2.5, 4); // 068 FL_BKstarmumu_2.5_4
        add_bin(O::A_FB_B__KSTAR_MU_MU, 2.5, 4); // 069 AFB_BKstarmumu_2.5_4
        add_bin(O::S_3_B__KSTAR_MU_MU, 2.5, 4); // 070 S3_BKstarmumu_2.5_4
        add_bin(O::S_4_B__KSTAR_MU_MU, 2.5, 4); // 071 S4_BKstarmumu_2.5_4
        add_bin(O::S_5_B__KSTAR_MU_MU, 2.5, 4); // 072 S5_BKstarmumu_2.5_4
        add_bin(O::S_7_B__KSTAR_MU_MU, 2.5, 4); // 073 S7_BKstarmumu_2.5_4
        add_bin(O::S_8_B__KSTAR_MU_MU, 2.5, 4); // 074 S8_BKstarmumu_2.5_4
        add_bin(O::S_9_B__KSTAR_MU_MU, 2.5, 4); // 075 S9_BKstarmumu_2.5_4
        add_bin(O::F_L_B__KSTAR_MU_MU, 4, 6); // 076 FL_BKstarmumu_4_6
        add_bin(O::A_FB_B__KSTAR_MU_MU, 4, 6); // 077 AFB_BKstarmumu_4_6
        add_bin(O::S_3_B__KSTAR_MU_MU, 4, 6); // 078 S3_BKstarmumu_4_6
        add_bin(O::S_4_B__KSTAR_MU_MU, 4, 6); // 079 S4_BKstarmumu_4_6
        add_bin(O::S_5_B__KSTAR_MU_MU, 4, 6); // 080 S5_BKstarmumu_4_6
        add_bin(O::S_7_B__KSTAR_MU_MU, 4, 6); // 081 S7_BKstarmumu_4_6
        add_bin(O::S_8_B__KSTAR_MU_MU, 4, 6); // 082 S8_BKstarmumu_4_6
        add_bin(O::S_9_B__KSTAR_MU_MU, 4, 6); // 083 S9_BKstarmumu_4_6
        add_bin(O::F_L_B__KSTAR_MU_MU, 15, 17); // 084 FL_BKstarmumu_15_17
        add_bin(O::A_FB_B__KSTAR_MU_MU, 15, 17); // 085 AFB_BKstarmumu_15_17
        add_bin(O::S_3_B__KSTAR_MU_MU, 15, 17); // 086 S3_BKstarmumu_15_17
        add_bin(O::S_4_B__KSTAR_MU_MU, 15, 17); // 087 S4_BKstarmumu_15_17
        add_bin(O::S_5_B__KSTAR_MU_MU, 15, 17); // 088 S5_BKstarmumu_15_17
        add_bin(O::S_7_B__KSTAR_MU_MU, 15, 17); // 089 S7_BKstarmumu_15_17
        add_bin(O::S_8_B__KSTAR_MU_MU, 15, 17); // 090 S8_BKstarmumu_15_17
        add_bin(O::S_9_B__KSTAR_MU_MU, 15, 17); // 091 S9_BKstarmumu_15_17
        add_bin(O::F_L_B__KSTAR_MU_MU, 17, 19); // 092 FL_BKstarmumu_17_19
        add_bin(O::A_FB_B__KSTAR_MU_MU, 17, 19); // 093 AFB_BKstarmumu_17_19
        add_bin(O::S_3_B__KSTAR_MU_MU, 17, 19); // 094 S3_BKstarmumu_17_19
        add_bin(O::S_4_B__KSTAR_MU_MU, 17, 19); // 095 S4_BKstarmumu_17_19
        add_bin(O::S_5_B__KSTAR_MU_MU, 17, 19); // 096 S5_BKstarmumu_17_19
        add_bin(O::S_7_B__KSTAR_MU_MU, 17, 19); // 097 S7_BKstarmumu_17_19
        add_bin(O::S_8_B__KSTAR_MU_MU, 17, 19); // 098 S8_BKstarmumu_17_19
        add_bin(O::S_9_B__KSTAR_MU_MU, 17, 19); // 099 S9_BKstarmumu_17_19
        add_bin(O::R_1_B__KSTAR_L_L, 0.045, 6); // 100 R-1_BKstarll_0.045_6
        // // TODO missing enum/mapping: 101 R-1_B0K0ll_1.1_6
        add_bin(O::R_1_B0__K0_L_L, 1.1, 6);
        add_bin(O::R_1_B__K_L_L, 1, 6); // 102 R-1_BKll_1_6_Belle
        add_bin(O::F_H_B__K_MU_MU, 1, 6); // 103 FH_BKmumu_1_6_CMS


        add_bin(O::DBR_DQ2_B0__KSTAR0_E_E, 0.0009, 1); // 104 dGamma/dq2_B0Kstar0ee_0.0009_1 //[2026-05-24_01] [WARN] Rejected MC nuisance sample 1 while trying to fill accepted sample 1 of 1000 : MC prediction contains non-finite observable
        add_bin(O::F_L_B0__KSTAR0_E_E, 0.0008, 0.257); // 105 FL_B0Kstar0ee_0.0008_0.257 //Rejected MC nuisance sample 125 while trying to fill accepted sample 1 of 1000 : MC prediction contains non-finite observable
        add_bin(O::A_T_RE_B0__KSTAR0_E_E, 0.0008, 0.257); // 106 ATRe_B0Kstar0ee_0.0008_0.257 //[2026-05-24_01] [WARN] Rejected MC nuisance sample 1 while trying to fill accepted sample 1 of 1000 : MC prediction contains non-finite observable
        add_bin(O::A_T_2_B0__KSTAR0_E_E, 0.0008, 0.257); // 107 AT2_B0Kstar0ee_0.0008_0.257 // [2026-05-24_01] [WARN] Rejected MC nuisance sample 1 while trying to fill accepted sample 1 of 1000 : MC prediction contains non-finite observable
        add_bin(O::DBR_DQ2_B__K_MU_MU, 0.1, 0.98); // 108 dGamma/dq2_BKmumu_0.1_0.98_CMS
        add_bin(O::DBR_DQ2_B__K_MU_MU, 1.1, 2); // 109 dGamma/dq2_BKmumu_1.1_2_CMS
        add_bin(O::DBR_DQ2_B__K_MU_MU, 2, 3); // 110 dGamma/dq2_BKmumu_2_3_CMS
        add_bin(O::DBR_DQ2_B__K_MU_MU, 3, 4); // 111 dGamma/dq2_BKmumu_3_4_CMS
        add_bin(O::DBR_DQ2_B__K_MU_MU, 4, 5); // 112 dGamma/dq2_BKmumu_4_5_CMS
        add_bin(O::DBR_DQ2_B__K_MU_MU, 5, 6); // 113 dGamma/dq2_BKmumu_5_6_CMS
        add_bin(O::DBR_DQ2_B__K_MU_MU, 14.82, 16); // 114 dGamma/dq2_BKmumu_14.82_16_CMS
        add_bin(O::DBR_DQ2_B__K_MU_MU, 16, 17); // 115 dGamma/dq2_BKmumu_16_17_CMS
        add_bin(O::DBR_DQ2_B__K_MU_MU, 17, 18); // 116 dGamma/dq2_BKmumu_17_18_CMS
        add_bin(O::DBR_DQ2_B__K_MU_MU, 18, 19.24); // 117 dGamma/dq2_BKmumu_18_19.24_CMS
        add_bin(O::DBR_DQ2_B__K_MU_MU, 19.24, 22.9); // 118 dGamma/dq2_BKmumu_19.24_22.9_CMS

        //FIRST STOP
        // add_bin(O::R_1_B__K_L_L, 1.1, 6); // 119 R-1_BKll_1.1_6_CMS
        add_bin(O::F_L_B0__KSTAR0_MU_MU, 1.1, 2); // 120 FL_B0Kstar0mumu_1.1_2_CMS
        add_bin(O::P_1_B0__KSTAR0_MU_MU, 1.1, 2); // 121 P1_B0Kstar0mumu_1.1_2_CMS
        add_bin(O::P_2_B0__KSTAR0_MU_MU, 1.1, 2); // 122 P2_B0Kstar0mumu_1.1_2_CMS
        add_bin(O::P_3_B0__KSTAR0_MU_MU, 1.1, 2); // 123 P3_B0Kstar0mumu_1.1_2_CMS
        add_bin(O::P_PRIME_4_B0__KSTAR0_MU_MU, 1.1, 2); // 124 P4prime_B0Kstar0mumu_1.1_2_CMS
        add_bin(O::P_PRIME_5_B0__KSTAR0_MU_MU, 1.1, 2); // 125 P5prime_B0Kstar0mumu_1.1_2_CMS
        add_bin(O::P_PRIME_6_B0__KSTAR0_MU_MU, 1.1, 2); // 126 P6prime_B0Kstar0mumu_1.1_2_CMS
        add_bin(O::P_PRIME_8_B0__KSTAR0_MU_MU, 1.1, 2); // 127 P8prime_B0Kstar0mumu_1.1_2_CMS
        add_bin(O::F_L_B0__KSTAR0_MU_MU, 2, 4.3); // 128 FL_B0Kstar0mumu_2_4.3_CMS
        add_bin(O::P_1_B0__KSTAR0_MU_MU, 2, 4.3); // 129 P1_B0Kstar0mumu_2_4.3_CMS
        add_bin(O::P_2_B0__KSTAR0_MU_MU, 2, 4.3); // 130 P2_B0Kstar0mumu_2_4.3_CMS
        add_bin(O::P_3_B0__KSTAR0_MU_MU, 2, 4.3); // 131 P3_B0Kstar0mumu_2_4.3_CMS
        add_bin(O::P_PRIME_4_B0__KSTAR0_MU_MU, 2, 4.3); // 132 P4prime_B0Kstar0mumu_2_4.3_CMS
        add_bin(O::P_PRIME_5_B0__KSTAR0_MU_MU, 2, 4.3); // 133 P5prime_B0Kstar0mumu_2_4.3_CMS
        add_bin(O::P_PRIME_6_B0__KSTAR0_MU_MU, 2, 4.3); // 134 P6prime_B0Kstar0mumu_2_4.3_CMS
        add_bin(O::P_PRIME_8_B0__KSTAR0_MU_MU, 2, 4.3); // 135 P8prime_B0Kstar0mumu_2_4.3_CMS
        add_bin(O::F_L_B0__KSTAR0_MU_MU, 4.3, 6); // 136 FL_B0Kstar0mumu_4.3_6_CMS
        add_bin(O::P_1_B0__KSTAR0_MU_MU, 4.3, 6); // 137 P1_B0Kstar0mumu_4.3_6_CMS
        add_bin(O::P_2_B0__KSTAR0_MU_MU, 4.3, 6); // 138 P2_B0Kstar0mumu_4.3_6_CMS
        add_bin(O::P_3_B0__KSTAR0_MU_MU, 4.3, 6); // 139 P3_B0Kstar0mumu_4.3_6_CMS
        add_bin(O::P_PRIME_4_B0__KSTAR0_MU_MU, 4.3, 6); // 140 P4prime_B0Kstar0mumu_4.3_6_CMS
        add_bin(O::P_PRIME_5_B0__KSTAR0_MU_MU, 4.3, 6); // 141 P5prime_B0Kstar0mumu_4.3_6_CMS
        add_bin(O::P_PRIME_6_B0__KSTAR0_MU_MU, 4.3, 6); // 142 P6prime_B0Kstar0mumu_4.3_6_CMS
        add_bin(O::P_PRIME_8_B0__KSTAR0_MU_MU, 4.3, 6); // 143 P8prime_B0Kstar0mumu_4.3_6_CMS
        add_bin(O::F_L_B0__KSTAR0_MU_MU, 14.18, 16); // 144 FL_B0Kstar0mumu_14.18_16_CMS
        add_bin(O::P_1_B0__KSTAR0_MU_MU, 14.18, 16); // 145 P1_B0Kstar0mumu_14.18_16_CMS
        add_bin(O::P_2_B0__KSTAR0_MU_MU, 14.18, 16); // 146 P2_B0Kstar0mumu_14.18_16_CMS
        add_bin(O::P_3_B0__KSTAR0_MU_MU, 14.18, 16); // 147 P3_B0Kstar0mumu_14.18_16_CMS
        add_bin(O::P_PRIME_4_B0__KSTAR0_MU_MU, 14.18, 16); // 148 P4prime_B0Kstar0mumu_14.18_16_CMS
        add_bin(O::P_PRIME_5_B0__KSTAR0_MU_MU, 14.18, 16); // 149 P5prime_B0Kstar0mumu_14.18_16_CMS
        add_bin(O::P_PRIME_6_B0__KSTAR0_MU_MU, 14.18, 16); // 150 P6prime_B0Kstar0mumu_14.18_16_CMS
        add_bin(O::P_PRIME_8_B0__KSTAR0_MU_MU, 14.18, 16); // 151 P8prime_B0Kstar0mumu_14.18_16_CMS
        add_bin(O::DBR_DQ2_BS__PHI_E_E, 0.1, 1.1); // 152 dGamma/dq2_Bsphiee_0.1_1.1
        add_bin(O::DBR_DQ2_BS__PHI_E_E, 1.1, 6); // 153 dGamma/dq2_Bsphiee_1.1_6
        add_bin(O::DBR_DQ2_BS__PHI_E_E, 15, 19); // 154 dGamma/dq2_Bsphiee_15_19 //[2026-05-24_01] [WARN] Rejected MC nuisance sample 1 while trying to fill accepted sample 1 of 1000 : MC prediction contains non-finite observable
        add_bin(O::R_1_BS__PHI_L_L, 0.1, 1.1); // 155 R-1_Bsphill_0.1_1.1
        add_bin(O::R_1_BS__PHI_L_L, 1.1, 6); // 156 R-1_Bsphill_1.1_6
        add_bin(O::R_1_BS__PHI_L_L, 15, 19); // 157 R-1_Bsphill_15_19 //[2026-05-24_01] [WARN] Rejected MC nuisance sample 1 while trying to fill accepted sample 1 of 1000 : MC prediction contains non-finite observable
        add_bin(O::F_L_BS_PHI_E_E, 0.0009, 0.2615); // 158 FL_Bsphiee_0.0009_0.2615 //[2026-05-24_01] [WARN] Rejected MC nuisance sample 1 while trying to fill accepted sample 1 of 1000 : MC prediction contains non-finite observable
        add_bin(O::A_T_2_BS_PHI_E_E, 0.0009, 0.2615); // 159 AT2_Bsphiee_0.0009_0.2615 //[2026-05-24_01] [WARN] Rejected MC nuisance sample 1 while trying to fill accepted sample 1 of 1000 : MC prediction contains non-finite observable
        add_bin(O::A_T_2_B0__KSTAR0_E_E, 0.0008, 1.12); // 160 AT2_B0Kstar0ee_0.0008_1.12_Belle //[2026-05-24_01] [WARN] Rejected MC nuisance sample 1 while trying to fill accepted sample 1 of 1000 : MC prediction contains non-finite observable
        add_bin(O::F_L_B0__KSTAR0_E_E, 1.1, 6); // 161 FL_B0Kstar0ee_1.1_6
        add_bin(O::P_1_B0__KSTAR0_E_E, 1.1, 6); // 162 P1_B0Kstar0ee_1.1_6
        add_bin(O::P_2_B0__KSTAR0_E_E, 1.1, 6); // 163 P2_B0Kstar0ee_1.1_6
        add_bin(O::P_3_B0__KSTAR0_E_E, 1.1, 6); // 164 P3_B0Kstar0ee_1.1_6
        add_bin(O::P_PRIME_4_B0__KSTAR0_E_E, 1.1, 6); // 165 P4prime_B0Kstar0ee_1.1_6
        add_bin(O::P_PRIME_5_B0__KSTAR0_E_E, 1.1, 6); // 166 P5prime_B0Kstar0ee_1.1_6
        add_bin(O::P_PRIME_6_B0__KSTAR0_E_E, 1.1, 6); // 167 P6prime_B0Kstar0ee_1.1_6
        add_bin(O::P_PRIME_8_B0__KSTAR0_E_E, 1.1, 6); // 168 P8prime_B0Kstar0ee_1.1_6
        add_bin(O::F_L_B0__KSTAR0_MU_MU, 0.06, 0.98); // 169 FL_B0Kstar0mumu_0.06_0.98_LHCb2025c2
        add_bin(O::S_2S_B0__KSTAR0_MU_MU, 0.06, 0.98); // 170 S2s_B0Kstar0mumu_0.06_0.98_LHCb2025c2
        add_bin(O::S_1C_B0__KSTAR0_MU_MU, 0.06, 0.98); // 171 S1c_B0Kstar0mumu_0.06_0.98_LHCb2025c2
        add_bin(O::P_1_B0__KSTAR0_MU_MU, 0.06, 0.98); // 172 P1_B0Kstar0mumu_0.06_0.98_LHCb2025c2
        add_bin(O::P_2_B0__KSTAR0_MU_MU, 0.06, 0.98); // 173 P2_B0Kstar0mumu_0.06_0.98_LHCb2025c2
        add_bin(O::P_3_B0__KSTAR0_MU_MU, 0.06, 0.98); // 174 P3_B0Kstar0mumu_0.06_0.98_LHCb2025c2
        add_bin(O::P_PRIME_4_B0__KSTAR0_MU_MU, 0.06, 0.98); // 175 P4prime_B0Kstar0mumu_0.06_0.98_LHCb2025c2
        add_bin(O::P_PRIME_5_B0__KSTAR0_MU_MU, 0.06, 0.98); // 176 P5prime_B0Kstar0mumu_0.06_0.98_LHCb2025c2
        add_bin(O::P_PRIME_6_B0__KSTAR0_MU_MU, 0.06, 0.98); // 177 P6prime_B0Kstar0mumu_0.06_0.98_LHCb2025c2
        add_bin(O::P_PRIME_8_B0__KSTAR0_MU_MU, 0.06, 0.98); // 178 P8prime_B0Kstar0mumu_0.06_0.98_LHCb2025c2
        add_bin(O::S_6C_B0__KSTAR0_MU_MU, 0.06, 0.98); // 179 S6c_B0Kstar0mumu_0.06_0.98_LHCb2025c2
        add_bin(O::DBR_DQ2_B0__KSTAR0_MU_MU, 0.06, 0.98); // 180 dGamma/dq2_B0Kstar0mumu_0.06_0.98_LHCb2025c2
        add_bin(O::F_L_B0__KSTAR0_MU_MU, 1.1, 2.5); // 181 FL_B0Kstar0mumu_1.1_2.5_LHCb2025c2
        add_bin(O::S_1C_B0__KSTAR0_MU_MU, 1.1, 2.5); // 182 S1c_B0Kstar0mumu_1.1_2.5_LHCb2025c2
        add_bin(O::P_1_B0__KSTAR0_MU_MU, 1.1, 2.5); // 183 P1_B0Kstar0mumu_1.1_2.5_LHCb2025c2
        add_bin(O::P_2_B0__KSTAR0_MU_MU, 1.1, 2.5); // 184 P2_B0Kstar0mumu_1.1_2.5_LHCb2025c2
        add_bin(O::P_3_B0__KSTAR0_MU_MU, 1.1, 2.5); // 185 P3_B0Kstar0mumu_1.1_2.5_LHCb2025c2
        add_bin(O::P_PRIME_4_B0__KSTAR0_MU_MU, 1.1, 2.5); // 186 P4prime_B0Kstar0mumu_1.1_2.5_LHCb2025c2
        add_bin(O::P_PRIME_5_B0__KSTAR0_MU_MU, 1.1, 2.5); // 187 P5prime_B0Kstar0mumu_1.1_2.5_LHCb2025c2
        add_bin(O::P_PRIME_6_B0__KSTAR0_MU_MU, 1.1, 2.5); // 188 P6prime_B0Kstar0mumu_1.1_2.5_LHCb2025c2
        add_bin(O::P_PRIME_8_B0__KSTAR0_MU_MU, 1.1, 2.5); // 189 P8prime_B0Kstar0mumu_1.1_2.5_LHCb2025c2
        add_bin(O::DBR_DQ2_B0__KSTAR0_MU_MU, 1.1, 2.5); // 190 dGamma/dq2_B0Kstar0mumu_1.1_2.5_LHCb2025c2
        add_bin(O::F_L_B0__KSTAR0_MU_MU, 2.5, 4); // 191 FL_B0Kstar0mumu_2.5_4.0_LHCb2025c2
        add_bin(O::S_1C_B0__KSTAR0_MU_MU, 2.5, 4); // 192 S1c_B0Kstar0mumu_2.5_4.0_LHCb2025c2
        add_bin(O::P_1_B0__KSTAR0_MU_MU, 2.5, 4); // 193 P1_B0Kstar0mumu_2.5_4.0_LHCb2025c2
        add_bin(O::P_2_B0__KSTAR0_MU_MU, 2.5, 4); // 194 P2_B0Kstar0mumu_2.5_4.0_LHCb2025c2
        add_bin(O::P_3_B0__KSTAR0_MU_MU, 2.5, 4); // 195 P3_B0Kstar0mumu_2.5_4.0_LHCb2025c2
        add_bin(O::P_PRIME_4_B0__KSTAR0_MU_MU, 2.5, 4); // 196 P4prime_B0Kstar0mumu_2.5_4.0_LHCb2025c2
        add_bin(O::P_PRIME_5_B0__KSTAR0_MU_MU, 2.5, 4); // 197 P5prime_B0Kstar0mumu_2.5_4.0_LHCb2025c2
        add_bin(O::P_PRIME_6_B0__KSTAR0_MU_MU, 2.5, 4); // 198 P6prime_B0Kstar0mumu_2.5_4.0_LHCb2025c2
        add_bin(O::P_PRIME_8_B0__KSTAR0_MU_MU, 2.5, 4); // 199 P8prime_B0Kstar0mumu_2.5_4.0_LHCb2025c2
        add_bin(O::DBR_DQ2_B0__KSTAR0_MU_MU, 2.5, 4); // 200 dGamma/dq2_B0Kstar0mumu_2.5_4.0_LHCb2025c2
        add_bin(O::F_L_B0__KSTAR0_MU_MU, 4, 6); // 201 FL_B0Kstar0mumu_4.0_6.0_LHCb2025c2
        add_bin(O::S_1C_B0__KSTAR0_MU_MU, 4, 6); // 202 S1c_B0Kstar0mumu_4.0_6.0_LHCb2025c2
        add_bin(O::P_1_B0__KSTAR0_MU_MU, 4, 6); // 203 P1_B0Kstar0mumu_4.0_6.0_LHCb2025c2
        add_bin(O::P_2_B0__KSTAR0_MU_MU, 4, 6); // 204 P2_B0Kstar0mumu_4.0_6.0_LHCb2025c2
        add_bin(O::P_3_B0__KSTAR0_MU_MU, 4, 6); // 205 P3_B0Kstar0mumu_4.0_6.0_LHCb2025c2
        add_bin(O::P_PRIME_4_B0__KSTAR0_MU_MU, 4, 6); // 206 P4prime_B0Kstar0mumu_4.0_6.0_LHCb2025c2
        add_bin(O::P_PRIME_5_B0__KSTAR0_MU_MU, 4, 6); // 207 P5prime_B0Kstar0mumu_4.0_6.0_LHCb2025c2
        add_bin(O::P_PRIME_6_B0__KSTAR0_MU_MU, 4, 6); // 208 P6prime_B0Kstar0mumu_4.0_6.0_LHCb2025c2
        add_bin(O::P_PRIME_8_B0__KSTAR0_MU_MU, 4, 6); // 209 P8prime_B0Kstar0mumu_4.0_6.0_LHCb2025c2
        add_bin(O::DBR_DQ2_B0__KSTAR0_MU_MU, 4, 6); // 210 dGamma/dq2_B0Kstar0mumu_4.0_6.0_LHCb2025c2
        add_bin(O::F_L_B0__KSTAR0_MU_MU, 15, 17); // 211 FL_B0Kstar0mumu_15.0_17.0_LHCb2025c2
        add_bin(O::S_1C_B0__KSTAR0_MU_MU, 15, 17); // 212 S1c_B0Kstar0mumu_15.0_17.0_LHCb2025c2
        add_bin(O::P_1_B0__KSTAR0_MU_MU, 15, 17); // 213 P1_B0Kstar0mumu_15.0_17.0_LHCb2025c2
        add_bin(O::P_2_B0__KSTAR0_MU_MU, 15, 17); // 214 P2_B0Kstar0mumu_15.0_17.0_LHCb2025c2
        add_bin(O::P_3_B0__KSTAR0_MU_MU, 15, 17); // 215 P3_B0Kstar0mumu_15.0_17.0_LHCb2025c2
        add_bin(O::P_PRIME_4_B0__KSTAR0_MU_MU, 15, 17); // 216 P4prime_B0Kstar0mumu_15.0_17.0_LHCb2025c2
        add_bin(O::P_PRIME_5_B0__KSTAR0_MU_MU, 15, 17); // 217 P5prime_B0Kstar0mumu_15.0_17.0_LHCb2025c2
        add_bin(O::P_PRIME_6_B0__KSTAR0_MU_MU, 15, 17); // 218 P6prime_B0Kstar0mumu_15.0_17.0_LHCb2025c2
        add_bin(O::P_PRIME_8_B0__KSTAR0_MU_MU, 15, 17); // 219 P8prime_B0Kstar0mumu_15.0_17.0_LHCb2025c2
        add_bin(O::DBR_DQ2_B0__KSTAR0_MU_MU, 15, 17); // 220 dGamma/dq2_B0Kstar0mumu_15.0_17.0_LHCb2025c2
        add_bin(O::F_L_B0__KSTAR0_MU_MU, 17, 19); // 221 FL_B0Kstar0mumu_17.0_19.0_LHCb2025c2
        add_bin(O::S_1C_B0__KSTAR0_MU_MU, 17, 19); // 222 S1c_B0Kstar0mumu_17.0_19.0_LHCb2025c2
        add_bin(O::P_1_B0__KSTAR0_MU_MU, 17, 19); // 223 P1_B0Kstar0mumu_17.0_19.0_LHCb2025c2
        add_bin(O::P_2_B0__KSTAR0_MU_MU, 17, 19); // 224 P2_B0Kstar0mumu_17.0_19.0_LHCb2025c2
        add_bin(O::P_3_B0__KSTAR0_MU_MU, 17, 19); // 225 P3_B0Kstar0mumu_17.0_19.0_LHCb2025c2
        add_bin(O::P_PRIME_4_B0__KSTAR0_MU_MU, 17, 19); // 226 P4prime_B0Kstar0mumu_17.0_19.0_LHCb2025c2
        add_bin(O::P_PRIME_5_B0__KSTAR0_MU_MU, 17, 19); // 227 P5prime_B0Kstar0mumu_17.0_19.0_LHCb2025c2
        add_bin(O::P_PRIME_6_B0__KSTAR0_MU_MU, 17, 19); // 228 P6prime_B0Kstar0mumu_17.0_19.0_LHCb2025c2
        add_bin(O::P_PRIME_8_B0__KSTAR0_MU_MU, 17, 19); // 229 P8prime_B0Kstar0mumu_17.0_19.0_LHCb2025c2
        add_bin(O::DBR_DQ2_B0__KSTAR0_MU_MU, 17, 19); // 230 dGamma/dq2_B0Kstar0mumu_17.0_19.0_LHCb2025c2

    std::shared_ptr<IStatParamOptimizerProxy> spop = std::make_shared<StatParamOptimizerProxy>();
    auto model = std::make_shared<ObservableInterfaceProxy>(oint, spop);
    
    // LOG_INFO(oint->compute_observable(Observables::IA_B__KSTAR_GAMMA)[0].value);
    // LOG_INFO(oint->get_exp_value(Observables::IA_B__KSTAR_GAMMA));
    // LOG_INFO(oint->compute_observable(Observables::BR_B_XS_GAMMA)[0].value);
    // LOG_INFO(oint->compute_observable(Observables::BR_BS_EE_UNTAG)[0].value);
    // LOG_INFO(oint->compute_observable(Observables::BR_BS_MUMU_UNTAG)[0].value);
    // LOG_INFO(oint->compute_observable(Observables::BR_B__Xs_mu_mu)[0].value);
    // LOG_INFO(oint->compute_observable(Observables::BR_B__Xs_mu_mu)[1].value);
    // LOG_INFO(oint->compute_observable(Observables::BR_B__Xs_e_e)[0].value);
    // LOG_INFO(oint->compute_observable(Observables::BR_B__Xs_e_e)[1].value);
    // LOG_INFO(oint->compute_observable(Observables::BR_B0__KSTAR0_GAMMA)[0].value);
    // LOG_INFO(oint->compute_observable(Observables::BR_B__KSTAR_GAMMA)[0].value);
    // LOG_INFO(oint->compute_observable(Observables::DBR_DQ2_B__KSTAR_MU_MU)[0].value);
    // LOG_INFO(oint->compute_observable(Observables::DBR_DQ2_B__KSTAR_MU_MU)[1].value);

    // exit(0);

    StatisticConfig config;
    config.advanced.chi2_covariance_ridge_rel = 1e-8;
    config.advanced.chi2_covariance_ridge_abs = 0.0;
    config.advanced.MLE_max_iter = 120000;
    config.advanced.MLE_tol = 0.01;
    config.advanced.MLE_trace_first_evals  = true;
    config.advanced.MLE_trace_max_evals  = 20;
    config.advanced.likelihood_mode = StatisticLikelihoodMode::CHI2_MC_COVARIANCE;
    config.MC_draws = 100;
    config.advanced.nuisance_sensitivity_contexts = -1;
    const std::string had_bsm_block =
        GroupMapper::str(WGroup::B, ScaleType::HADRONIC, WilsonBasis::B_STANDARD)
        + "__BSM_INTERMEDIATE";

    const std::string had_bsm_block2 =
        GroupMapper::str(WGroup::BScalar, ScaleType::HADRONIC, WilsonBasis::B_STANDARD)
        + "__BSM_INTERMEDIATE";

    std::vector<ParamId> p_specs = {
        ParamId{ParameterType::WILSON, had_bsm_block, WCoefMapper::flha_full(WCoef::C9, QCDOrder::LO, ContributionType::BSM)},
        ParamId{ParameterType::WILSON, had_bsm_block, WCoefMapper::flha_full(WCoef::C10, QCDOrder::LO, ContributionType::BSM)},
        // ParamId{ParameterType::WILSON, had_bsm_block, WCoefMapper::flha_full(WCoef::C7, QCDOrder::LO, ContributionType::BSM)},
        // ParamId{ParameterType::WILSON, had_bsm_block, WCoefMapper::flha_full(WCoef::C8, QCDOrder::LO, ContributionType::BSM)},
        // ParamId{ParameterType::WILSON, had_bsm_block2, WCoefMapper::flha_full(WCoef::CQ1, QCDOrder::LO, ContributionType::BSM)},
        // ParamId{ParameterType::WILSON, had_bsm_block2, WCoefMapper::flha_full(WCoef::CQ2, QCDOrder::LO, ContributionType::BSM)},
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
    // std::set<std::string> exp = {"CMS", "LHCb2025c2"};
    std::set<std::string> exp = {"DEFAULT", "Belle", "CMS", "LHCb2025c2"};
    // std::set<std::string> exp = {"LHCb2020"};
    // std::set<std::string> exp = {"LHCb2025c2"};
    // std::set<std::string> exp = {"CMS"};
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

        // std::array<double, 4> bounds = {-20, 20, -20, 20};
        std::array<double, 4> bounds = {
            -3.0, 2.0,
            -2.5, 2.5
        };

        auto trace = std::make_shared<std::ofstream>("contour_trace.csv");
        (*trace) << "type,level,path_id,point_id,x,y,n_paths,n_points,elapsed_s,message\n";

        ContourOptions  opt;
        opt.primary_contour_method = ContourAlgorithm::MINUIT;
        opt.fallback_contour_method = ContourAlgorithm::AMS;
        opt.profile_backend = ProfileBackend::MINUIT;
        opt.resolution = 10;

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
        auto t4 = std::chrono::steady_clock::now();
        auto c95 = stat.confidence_contour(p1, p2, 2, bounds, opt);
        auto t5 = std::chrono::steady_clock::now();
        
        auto c68 = stat.confidence_contour(p1, p2, 1, bounds, opt);
        std::cout << "finish first contour" << std::endl;


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
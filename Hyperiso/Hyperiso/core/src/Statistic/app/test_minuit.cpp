#include "Logger.h"
#include <iostream>
#include <cassert>
#include "ObservableInterface.h"
#include "HyperisoMaster.h"
#include "config.hpp"
#include "BlockProxy.h"
#include "StatisticInterface.h"

struct ObservableUncertainty {
    double err_down; // erreur vers le bas
    double err_up;   // erreur vers le haut
    double err_sym; // erreur vers tout
};

static ObservableUncertainty uncertainty_from_summary(const GaussianSummary& gs) {
    ObservableUncertainty out{};

    if (gs.symmetric) {
        out.err_down = std::abs(gs.sigma);
        out.err_up   = std::abs(gs.sigma);
    } else {
        // Convention naturelle :
        // sigma_m = erreur vers le bas
        // sigma_p = erreur vers le haut
        out.err_down = std::abs(gs.sigma_m);
        out.err_up   = std::abs(gs.sigma_p);

        // Sécurité si jamais sigma_m/p sont nuls mais sigma existe
        if (out.err_down == 0.0 && out.err_up == 0.0 && gs.sigma != 0.0) {
            out.err_down = std::abs(gs.sigma);
            out.err_up   = std::abs(gs.sigma);
        }
    }
    out.err_sym = gs.sigma;
    return out;
}


static void write_observables_to_csv(
    const std::string& filename,
    const std::vector<ObservableValue>& obs_val,
    const std::vector<ObservableUncertainty>& uncertainties
) {
    if (obs_val.size() != uncertainties.size()) {
        throw std::runtime_error("obs_val et uncertainties n'ont pas la même taille");
    }

    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Impossible d'ouvrir le fichier CSV: " + filename);
    }

    out << "bin_low,bin_high,value,err_down,err_up\n";
    out << std::setprecision(17);

    for (std::size_t i = 0; i < obs_val.size(); ++i) {
        const auto& obs = obs_val[i];
        const auto& unc = uncertainties[i];

        if (!obs.bin.has_value()) {
            continue;
        }

        out << obs.bin->first  << ","
            << obs.bin->second << ","
            << obs.value       << ","
            << unc.err_down    << ","
            << unc.err_up      << ","
            << unc.err_sym      << "\n";
    }
}

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);
    HyperisoMaster hyp;
    HyperisoConfig config;
    config.model = Model::SM;
    hyp.init("lha/si_input.flha", config);

    QCDOrder order = QCDOrder::NNLO;
    std::shared_ptr<ObservableInterface> oi = std::make_shared<ObservableInterface>();

    Decays dec = Decays::K__pi_nu_nu;
    // KllDecayConfig dec_cfg;
    // dec_cfg.gen = 1;
    // dec_cfg.N_L_sign = -1;

    // oi.set_decay_config(dec, dec_cfg);
    // oi.add_observable(Observables::TEST, order);
    // oi.compute_observable(Observables::TEST);
    std::vector<double> squares;

    for (double x = 0.05; x<8.1; x+=0.01) {
        squares.push_back(x);
    }
    for (double x = 15; x<20; x+=0.01) {
        squares.push_back(x);
    }
    for (auto elem : squares) {
        oi->add_observable(BinnedObservableId(ObservableMapper::to_id(Observables::F_L_B__KSTAR_MU_MU), {elem, elem+0.1}), QCDOrder::NNLO, false);
    }
    
    std::vector<ObservableValue> obs_val = oi->compute_observable(Observables::F_L_B__KSTAR_MU_MU);

    StatisticConfig sc = StatisticConfig();
    sc.MC_draws = 2000;
    StatisticInterface si = StatisticInterface(sc, oi);

    std::map<BinnedObservableId, GaussianSummary> unc_map = si.compute_uncertainties();

    std::vector<ObservableUncertainty> uncertainties;
    uncertainties.reserve(obs_val.size());

    for (const auto& obs : obs_val) {
        if (!obs.bin.has_value()) {
            uncertainties.push_back({0.0, 0.0});
            continue;
        }

        BinnedObservableId key(obs.id, *obs.bin);

        auto it = unc_map.find(key);
        std::cout
          << "OBS  [" << obs.bin->first << ", " << obs.bin->second << "] "
          << " value=" << obs.value
          << " | GS.id=[" << it->second.id.p.first << ", " << it->second.id.p.second << "]"
          << " mu=" << it->second.mu
          << " sigma=" << it->second.sigma
          << " sigma_m=" << it->second.sigma_m
          << " sigma_p=" << it->second.sigma_p
          << " skew=" << it->second.skew
          << " symmetric=" << it->second.symmetric
          << "\n";
        if (it == unc_map.end()) {
            std::cerr << "Warning: uncertainty not found for bin ["
                      << obs.bin->first << ", " << obs.bin->second << "]\n";

            // Tu peux choisir 0, ou bien throw si tu veux que ça casse tout de suite
            uncertainties.push_back({0.0, 0.0});
            continue;
        }

        uncertainties.push_back(uncertainty_from_summary(it->second));
    }

    write_observables_to_csv("K__pi_nu_nu.csv", obs_val, uncertainties);

    std::cout << "CSV écrit dans K__pi_nu_nu.csv\n";
    return 0;

    // oi.add_observables(dec, order, false);
    // for (auto o : DecayMapper::get_observables(dec)) {
    //     if (o == Observables::TEST) continue;

    //     // if (o == Observables::A_FB_B__KSTAR_L_L || o == Observables::F_L_B__KSTAR_L_L) {
    //         auto obs_values = oi.compute_observable(o);
    //     std::stringstream ss;
    //     ss << std::scientific << std::setprecision(3);
    //     if (obs_values.size() == 1) {
    //         ss << "= " << obs_values[0].value;
    //     } else {
    //         ss << ": ";
    //         for (auto ov : obs_values) {
    //             ss << "[" << ov.bin.value().first << ", " << ov.bin.value().second << "] = " << ov.value << ", ";
    //         }
    //     }
            
    //     LOG_INFO(ObservableMapper::str(o), ss.str());   
    //     // }
    // }

    // auto Gamma = oi.compute_observable(Observables::DGAMMA_DQ2_BS__PHI_L_L)[0].value;
    // auto Gamma_bar = oi.compute_observable(Observables::DGAMMA_BAR_DQ2_BS__PHI_L_L)[0].value;

    // std::stringstream ss;
    // ss << std::scientific << std::setprecision(3);
    // ss << "= " << (Gamma + Gamma_bar) / 2;
    // LOG_INFO("BR(Bs > phi mu mu)", ss.str());   

    return 0;
}
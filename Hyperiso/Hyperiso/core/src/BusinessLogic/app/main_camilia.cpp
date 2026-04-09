#include <chrono>
#include <iostream>

struct ScopedTimer {
    std::string name;
    std::chrono::steady_clock::time_point t0;
    ScopedTimer(std::string n) : name(std::move(n)), t0(std::chrono::steady_clock::now()) {}
    ~ScopedTimer() {
        auto t1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::cerr << "[TIMER] " << name << " : " << ms << " ms\n";
    }
};

#include "Logger.h"
#include <iostream>
#include <cassert>
#include "ObservableInterface.h"
#include "HyperisoMaster.h"
#include "config.hpp"
#include "BlockProxy.h"

struct ObservableUncertainty {
    double err_down; // erreur vers le bas
    double err_up;   // erreur vers le haut
};

void write_observables_to_csv(
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
            // Si jamais tu as des observables non binnées, on les skip ici.
            // Tu peux aussi décider d'écrire NaN à la place.
            continue;
        }

        const double bin_low  = obs.bin->first;
        const double bin_high = obs.bin->second;

        out << bin_low  << ","
            << bin_high << ","
            << obs.value << ","
            << unc.err_down << ","
            << unc.err_up   << "\n";
    }
}

int main() {
    ScopedTimer total("main total");
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);
    HyperisoMaster hyp;
    HyperisoConfig config;
    config.model = Model::SM;
    {
        ScopedTimer t("hyp.init");
        hyp.init("lha/si_input.flha", config);}

    QCDOrder order = QCDOrder::NNLO;
    ObservableInterface oi;

    Decays dec = Decays::K__pi_nu_nu;
    // KllDecayConfig dec_cfg;
    // dec_cfg.gen = 1;
    // dec_cfg.N_L_sign = -1;

    // oi.set_decay_config(dec, dec_cfg);
    // oi.add_observable(Observables::TEST, order);
    // oi.compute_observable(Observables::TEST);
    std::vector<double> squares;

    for (double x = 0.05; x<8.1; x+=1) {
        squares.push_back(x);
    }
    {
        ScopedTimer t("add_observable loop");
        for (auto elem : squares) {
        oi.add_observable(BinnedObservableId(ObservableMapper::to_id(Observables::DBR_DQ2_B__K_MU_MU), {elem, elem+0.1}), QCDOrder::NNLO, false);
        }
    }
    std::vector<ObservableValue> obs_val;
    {
        ScopedTimer t("compute_observable(F_L)");
        obs_val = oi.compute_observable(Observables::DBR_DQ2_B__K_MU_MU);
    }

    std::vector<ObservableUncertainty> uncertainties;
    uncertainties.reserve(obs_val.size());
    
    for (const auto& obs : obs_val) {
        double sigma = 0.10 * std::abs(obs.value);
        uncertainties.push_back({sigma, sigma});
    }

    {
        ScopedTimer t("write csv");
        write_observables_to_csv("DBR_DQ2_B__K_MU_MU.csv", obs_val, uncertainties);
    }

    std::cout << "CSV écrit dans DBR_DQ2_B__K_MU_MU.csv\n";
    return 0;
}
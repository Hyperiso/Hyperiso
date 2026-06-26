#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Include.h"
#include "Logger.h"
#include "ObservableInterface.h"
#include "HyperisoMaster.h"
#include "config.hpp"
#include "ParamOptimizerAdapter.h"

struct ScopedTimer {
    std::string name;
    std::chrono::steady_clock::time_point t0;
    explicit ScopedTimer(std::string n) : name(std::move(n)), t0(std::chrono::steady_clock::now()) {}
    ~ScopedTimer() {
        auto t1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::cerr << "[TIMER] " << name << " : " << ms << " ms\n";
    }
};

namespace {

struct Row {
    std::string si_name;       // nom côté SuperIso, pour merge/join direct des CSV
    Observables hyp_obs;
    std::string hyp_enum;
    bool binned;
    double q2min;
    double q2max;
};

std::string csv_escape(const std::string& s) {
    std::string out = "\"";
    for (char c : s) {
        if (c == '"') out += "\"\"";
        else out += c;
    }
    out += "\"";
    return out;
}

void write_num(std::ofstream& out, double x) {
    if (std::isfinite(x)) out << std::setprecision(17) << x;
    else out << "nan";
}

bool same_bin(double a, double b) {
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a - b) <= 1e-10 * scale;
}

BinnedObservableId make_id(const Row& r) {
    if (r.binned) return BinnedObservableId{ObservableMapper::to_id(r.hyp_obs), {r.q2min, r.q2max}};
    return BinnedObservableId{ObservableMapper::to_id(r.hyp_obs)};
}

double find_value(const std::map<ObservableId, std::vector<ObservableValue>>& all,
                  const Row& r) {
    const auto id = make_id(r);
    auto it = all.find(id.s);
    if (it == all.end()) {
        throw std::runtime_error("observable non trouvee dans compute_all(): " + r.hyp_enum);
    }

    for (const auto& v : it->second) {
        if (!r.binned) {
            if (!v.bin.has_value()) return v.value;
        } else if (v.bin.has_value()
                   && same_bin(v.bin->first, r.q2min)
                   && same_bin(v.bin->second, r.q2max)) {
            return v.value;
        }
    }

    throw std::runtime_error("bin non trouve pour " + r.hyp_enum);
}

void add_row(std::vector<Row>& rows,
             const std::string& si_name,
             Observables obs,
             const std::string& hyp_enum,
             double q2min = std::numeric_limits<double>::quiet_NaN(),
             double q2max = std::numeric_limits<double>::quiet_NaN()) {
    const bool binned = std::isfinite(q2min) && std::isfinite(q2max);
    rows.push_back(Row{si_name, obs, hyp_enum, binned, q2min, q2max});
}

std::vector<Row> make_rows() {
    using O = Observables;
    std::vector<Row> rows;

    add_row(rows, "AI_BKstargamma",       O::IA_B__KSTAR_GAMMA,      "IA_B__KSTAR_GAMMA");
    add_row(rows, "BR_BXsgamma",          O::BR_B_XS_GAMMA,          "BR_B_XS_GAMMA");
    add_row(rows, "BR_BXsmumu_1_6",       O::BR_B__Xs_mu_mu,         "BR_B__Xs_mu_mu", 1.0, 6.0);
    add_row(rows, "BR_BXsmumu_14.2_22",   O::BR_B__Xs_mu_mu,         "BR_B__Xs_mu_mu", 14.2, 22.0);
    add_row(rows, "BR_B0Kstar0gamma",     O::BR_B0__KSTAR0_GAMMA,    "BR_B0__KSTAR0_GAMMA");
    add_row(rows, "BR_BKstargamma",       O::BR_B__KSTAR_GAMMA,      "BR_B__KSTAR_GAMMA");

    const std::vector<std::pair<double,double>> bins = {
        {0.06, 0.98}, {1.1, 2.5}, {2.5, 4.0}, {4.0, 6.0}, {15.0, 17.0}, {17.0, 19.0}
    };

    auto tag = [](double x) {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(1) << x;
        std::string s = ss.str();
        if (s.size() >= 2 && s.substr(s.size() - 2) == ".0") s.erase(s.size() - 2);
        if (std::abs(x - 0.06) < 1e-12) return std::string("0.06");
        if (std::abs(x - 0.98) < 1e-12) return std::string("0.98");
        return s;
    };

    for (const auto& [lo, hi] : bins) {
        const std::string suffix = tag(lo) + "_" + tag(hi) + "_LHCb2025c2";
        add_row(rows, "FL_B0Kstar0mumu_"       + suffix, O::F_L_B0__KSTAR0_MU_MU,        "F_L_B0__KSTAR0_MU_MU", lo, hi);
        if (lo == 0.06) {
            add_row(rows, "S2s_B0Kstar0mumu_"  + suffix, O::S_2S_B0__KSTAR0_MU_MU,       "S_2S_B0__KSTAR0_MU_MU", lo, hi);
        }
        add_row(rows, "S1c_B0Kstar0mumu_"      + suffix, O::S_1C_B0__KSTAR0_MU_MU,       "S_1C_B0__KSTAR0_MU_MU", lo, hi);
        add_row(rows, "P1_B0Kstar0mumu_"       + suffix, O::P_1_B0__KSTAR0_MU_MU,        "P_1_B0__KSTAR0_MU_MU", lo, hi);
        add_row(rows, "P2_B0Kstar0mumu_"       + suffix, O::P_2_B0__KSTAR0_MU_MU,        "P_2_B0__KSTAR0_MU_MU", lo, hi);
        add_row(rows, "P3_B0Kstar0mumu_"       + suffix, O::P_3_B0__KSTAR0_MU_MU,        "P_3_B0__KSTAR0_MU_MU", lo, hi);
        add_row(rows, "P4prime_B0Kstar0mumu_"  + suffix, O::P_PRIME_4_B0__KSTAR0_MU_MU, "P_PRIME_4_B0__KSTAR0_MU_MU", lo, hi);
        add_row(rows, "P5prime_B0Kstar0mumu_"  + suffix, O::P_PRIME_5_B0__KSTAR0_MU_MU, "P_PRIME_5_B0__KSTAR0_MU_MU", lo, hi);
        add_row(rows, "P6prime_B0Kstar0mumu_"  + suffix, O::P_PRIME_6_B0__KSTAR0_MU_MU, "P_PRIME_6_B0__KSTAR0_MU_MU", lo, hi);
        add_row(rows, "P8prime_B0Kstar0mumu_"  + suffix, O::P_PRIME_8_B0__KSTAR0_MU_MU, "P_PRIME_8_B0__KSTAR0_MU_MU", lo, hi);
        if (lo == 0.06) {
            add_row(rows, "S6c_B0Kstar0mumu_"  + suffix, O::S_6C_B0__KSTAR0_MU_MU,       "S_6C_B0__KSTAR0_MU_MU", lo, hi);
        }

        /* Important : SuperIso appelle ces lignes dGamma/dq2.
         * Ici on utilise volontairement DGAMMA_DQ2 et pas DBR_DQ2, pour vérifier
         * explicitement si ton main HyperIso utilisait accidentellement dBR/dq2.
         */
        add_row(rows, "dGamma/dq2_B0Kstar0mumu_" + suffix,
                O::DBR_DQ2_B0__KSTAR0_MU_MU,
                "DGAMMA_DQ2_B0__KSTAR0_MU_MU", lo, hi);
    }

    return rows;
}

void apply_bsm_point(double dC7, double dC8, double dC9, double dC10) {
    if (std::abs(dC7) + std::abs(dC8) + std::abs(dC9) + std::abs(dC10) <= 1e-15) return;

    const std::string had_bsm_block =
        GroupMapper::str(WGroup::B, ScaleType::HADRONIC, WilsonBasis::B_STANDARD)
        + "__BSM_INTERMEDIATE";

    ParamOptimizerAdapter spop{{ParameterType::SM, ParameterType::FLAVOR, ParameterType::DECAY, ParameterType::WILSON}};
    auto set_wc = [&](WCoef c, double v) {
        spop.set_value(had_bsm_block,
                       WCoefMapper::flha_full(c, QCDOrder::LO, ContributionType::BSM),
                       v);
    };

    set_wc(WCoef::C7,  dC7);
    set_wc(WCoef::C8,  dC8);
    set_wc(WCoef::C9,  dC9);
    set_wc(WCoef::C10, dC10);
    spop.commit();
}

} // namespace

int main(int argc, char** argv) {
    ScopedTimer total("main total");

    const std::string out_path = (argc >= 2) ? argv[1] : "hyperiso_sensitive.csv";
    const std::string lha_path = (argc >= 3) ? argv[2] : "lha/si_input.flha";
    const double dC7  = (argc >= 4) ? std::atof(argv[3]) : 0.0;
    const double dC8  = (argc >= 5) ? std::atof(argv[4]) : 0.0;
    const double dC9  = (argc >= 6) ? std::atof(argv[5]) : 0.0;
    const double dC10 = (argc >= 7) ? std::atof(argv[6]) : 0.0;

    try {
        Logger::getInstance()->setLevel(Logger::LogLevel::INFO);
        std::cerr << "[INFO] LHA = " << lha_path << "\n";

        HyperisoMaster hyp;
        HyperisoConfig config;
        config.model = Model::SM;

        {
            ScopedTimer t("hyp.init");
            hyp.init(lha_path, config);
        }

        apply_bsm_point(dC7, dC8, dC9, dC10);

        ObservableInterface oi;
        constexpr bool add_deps = false;

        BKstarllConfig cfg_BKs;
        cfg_BKs.ff_src = BV_FF_Src::GRvDV;
        oi.set_decay_config(Decays::B__Kstar_l_l, cfg_BKs);
        oi.set_bkstarll_threads(4);

        BKstarGammaConfig cfg_BKsgamma;
        cfg_BKsgamma.ff_src = BV_FF_Src::GRvDV;
        oi.set_decay_config(Decays::B__Kstar_gamma, cfg_BKsgamma);

        BsPhiConfig cfg_BsPhi;
        cfg_BsPhi.ff_src = BV_FF_Src::GRvDV;
        oi.set_decay_config(Decays::Bs__phi_l_l, cfg_BsPhi);
        oi.set_bsphi_threads(4);

        BKllConfig cfg_BK;
        cfg_BK.ff_src = BP_FF_Src::GKvD_SR_LAT;
        oi.set_decay_config(Decays::B__K_l_l, cfg_BK);
        oi.set_bkll_threads(4);

        const auto rows = make_rows();
        std::set<Observables> seen_unbinned;
        std::set<std::tuple<Observables,double,double>> seen_binned;

        {
            ScopedTimer t("add_observable loop");
            for (const auto& r : rows) {
                if (r.binned) {
                    auto key = std::make_tuple(r.hyp_obs, r.q2min, r.q2max);
                    if (seen_binned.insert(key).second) {
                        oi.add_observable(BinnedObservableId{ObservableMapper::to_id(r.hyp_obs), {r.q2min, r.q2max}},
                                          QCDOrder::NNLO,
                                          add_deps);
                    }
                } else {
                    if (seen_unbinned.insert(r.hyp_obs).second) {
                        oi.add_observable(ObservableMapper::to_id(r.hyp_obs), QCDOrder::NNLO, add_deps);
                    }
                }
            }
        }

        std::map<ObservableId, std::vector<ObservableValue>> all;
        {
            ScopedTimer t("compute_all");
            all = oi.compute_all();
        }

        std::ofstream out(out_path);
        if (!out) throw std::runtime_error("Impossible d'ouvrir le CSV: " + out_path);

        out << "code,si_name,hyper_enum,bin_low,bin_high,value,err_down,err_up,dC7,dC8,dC9,dC10,lha,status\n";
        out << std::setprecision(17);

        for (const auto& r : rows) {
            double value = std::numeric_limits<double>::quiet_NaN();
            std::string status = "ok";
            try {
                value = find_value(all, r);
            } catch (const std::exception& e) {
                status = e.what();
            }

            out << "hyperiso," << csv_escape(r.si_name) << "," << csv_escape(r.hyp_enum) << ",";
            if (r.binned) {
                write_num(out, r.q2min); out << ","; write_num(out, r.q2max);
            } else {
                out << "nan,nan";
            }
            out << ",";
            write_num(out, value);
            out << ",nan,nan,";
            write_num(out, dC7); out << ",";
            write_num(out, dC8); out << ",";
            write_num(out, dC9); out << ",";
            write_num(out, dC10); out << ",";
            out << csv_escape(lha_path) << ",";
            out << csv_escape(status) << "\n";
        }

        std::cout << "CSV ecrit dans " << out_path << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "[FATAL] " << e.what() << "\n";
        return 1;
    }
}

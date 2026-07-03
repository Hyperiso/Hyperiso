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

/*
 * Diagnostic HyperIso pour les 230 observables de test_contour_BKstarmumu / myobs.in.
 *
 * Usage:
 *   ./hyperiso_allobs_compare.x [out.csv] [lha_path] [dC7] [dC8] [dC9] [dC10]
 *
 * Exemple:
 *   ./hyperiso_allobs_compare.x hyperiso_allobs_SM.csv lha/si_input.flha 0 0 0 0
 *
 * Les noms si_name sont exactement les noms SuperIso. Le CSV peut donc être mergé
 * directement avec la sortie de superiso_allobs_compare.c via compare_allobs.py.
 */

struct ScopedTimer {
    std::string name;
    std::chrono::steady_clock::time_point t0;
    explicit ScopedTimer(std::string n) : name(std::move(n)), t0(std::chrono::steady_clock::now()) {}
    ~ScopedTimer() {
        auto t1 = std::chrono::steady_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::cerr << "[TIMER] " << name << " : " << ms << " ms\n";
    }
};

namespace {

struct Row {
    int idx;
    std::string si_name;
    std::string experiment;
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
        throw std::runtime_error("observable absente de compute_all(): " + r.hyp_enum);
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

    throw std::runtime_error("bin absent pour " + r.hyp_enum
                             + " [" + std::to_string(r.q2min)
                             + ", " + std::to_string(r.q2max) + "]");
}

void add_row(std::vector<Row>& rows,
             int idx,
             const std::string& si_name,
             const std::string& experiment,
             Observables obs,
             const std::string& hyp_enum,
             double q2min = std::numeric_limits<double>::quiet_NaN(),
             double q2max = std::numeric_limits<double>::quiet_NaN()) {
    const bool binned = std::isfinite(q2min) && std::isfinite(q2max);
    rows.push_back(Row{idx, si_name, experiment, obs, hyp_enum, binned, q2min, q2max});
}

std::vector<Row> make_rows() {
    using O = Observables;
    std::vector<Row> rows;
    rows.reserve(230);

    add_row(rows, 1, "AI_BKstargamma", "DEFAULT", O::IA_B__KSTAR_GAMMA, "IA_B__KSTAR_GAMMA");
    add_row(rows, 2, "BR_BXsgamma", "DEFAULT", O::BR_B_XS_GAMMA, "BR_B_XS_GAMMA");
    add_row(rows, 3, "BRuntag_Bsmumu", "DEFAULT", O::BR_BS_MUMU_UNTAG, "BR_BS_MUMU_UNTAG");
    add_row(rows, 4, "BRuntag_Bsee", "DEFAULT", O::BR_BS_EE_UNTAG, "BR_BS_EE_UNTAG");
    add_row(rows, 5, "BR_BXsmumu_1_6", "DEFAULT", O::BR_B__Xs_mu_mu, "BR_B__Xs_mu_mu", 1, 6);
    add_row(rows, 6, "BR_BXsmumu_14.2_22", "DEFAULT", O::BR_B__Xs_mu_mu, "BR_B__Xs_mu_mu", 14.2, 22);
    add_row(rows, 7, "BR_BXsee_1_6", "DEFAULT", O::BR_B__Xs_e_e, "BR_B__Xs_e_e", 1, 6);
    add_row(rows, 8, "BR_BXsee_14.2_22", "DEFAULT", O::BR_B__Xs_e_e, "BR_B__Xs_e_e", 14.2, 22);
    add_row(rows, 9, "BR_B0Kstar0gamma", "DEFAULT", O::BR_B0__KSTAR0_GAMMA, "BR_B0__KSTAR0_GAMMA");
    add_row(rows, 10, "BR_BKstargamma", "DEFAULT", O::BR_B__KSTAR_GAMMA, "BR_B__KSTAR_GAMMA");
    add_row(rows, 11, "dGamma/dq2_BKstarmumu_1.1_6", "DEFAULT", O::DBR_DQ2_B__KSTAR_MU_MU, "DBR_DQ2_B__KSTAR_MU_MU", 1.1, 6);
    add_row(rows, 12, "dGamma/dq2_BKstarmumu_15_19", "DEFAULT", O::DBR_DQ2_B__KSTAR_MU_MU, "DBR_DQ2_B__KSTAR_MU_MU", 15, 19);
    add_row(rows, 13, "R-1_B0Kstar0ll_0.1_1.1", "DEFAULT", O::R_1_B0__KSTAR0_L_L, "R_1_B0__KSTAR0_L_L", 0.1, 1.1);
    add_row(rows, 14, "R-1_B0Kstar0ll_1.1_6", "DEFAULT", O::R_1_B0__KSTAR0_L_L, "R_1_B0__KSTAR0_L_L", 1.1, 6);
    add_row(rows, 15, "R-1_B0Kstar0ll_0.045_1.1_Belle", "Belle", O::R_1_B0__KSTAR0_L_L, "R_1_B0__KSTAR0_L_L", 0.045, 1.1);
    add_row(rows, 16, "R-1_B0Kstar0ll_1.1_6_Belle", "Belle", O::R_1_B0__KSTAR0_L_L, "R_1_B0__KSTAR0_L_L", 1.1, 6);
    add_row(rows, 17, "R-1_B0Kstar0ll_15_19_Belle", "Belle", O::R_1_B0__KSTAR0_L_L, "R_1_B0__KSTAR0_L_L", 15, 19);
    add_row(rows, 18, "dGamma/dq2_B0K0mumu_1.1_6", "DEFAULT", O::DBR_DQ2_B0__K0_MU_MU, "DBR_DQ2_B0__K0_MU_MU", 1.1, 6);
    add_row(rows, 19, "dGamma/dq2_B0K0mumu_15_22", "DEFAULT", O::DBR_DQ2_B0__K0_MU_MU, "DBR_DQ2_B0__K0_MU_MU", 15, 22);
    add_row(rows, 20, "dGamma/dq2_BKmumu_1.1_6", "DEFAULT", O::DBR_DQ2_B__K_MU_MU, "DBR_DQ2_B__K_MU_MU", 1.1, 6);
    add_row(rows, 21, "FH_BKmumu_1.1_6", "DEFAULT", O::F_H_B__K_MU_MU, "F_H_B__K_MU_MU", 1.1, 6);
    add_row(rows, 22, "dGamma/dq2_BKmumu_15_22", "DEFAULT", O::DBR_DQ2_B__K_MU_MU, "DBR_DQ2_B__K_MU_MU", 15, 22);
    add_row(rows, 23, "FH_BKmumu_15_22", "DEFAULT", O::F_H_B__K_MU_MU, "F_H_B__K_MU_MU", 15, 22);
    add_row(rows, 24, "R-1_BKll_0.1_1.1", "DEFAULT", O::R_1_B__K_L_L, "R_1_B__K_L_L", 0.1, 1.1);
    add_row(rows, 25, "R-1_BKll_1.1_6", "DEFAULT", O::R_1_B__K_L_L, "R_1_B__K_L_L", 1.1, 6);
    add_row(rows, 26, "dGamma/dq2_Bsphimumu_0.1_0.98", "DEFAULT", O::DBR_DQ2_BS__PHI_MU_MU, "DBR_DQ2_BS__PHI_MU_MU", 0.1, 0.98);
    add_row(rows, 27, "FL_Bsphimumu_0.1_0.98", "DEFAULT", O::F_L_BS_PHI_MU_MU, "F_L_BS_PHI_MU_MU", 0.1, 0.98);
    add_row(rows, 28, "S3_Bsphimumu_0.1_0.98", "DEFAULT", O::S_3_BS_PHI_MU_MU, "S_3_BS_PHI_MU_MU", 0.1, 0.98);
    add_row(rows, 29, "S4_Bsphimumu_0.1_0.98", "DEFAULT", O::S_4_BS_PHI_MU_MU, "S_4_BS_PHI_MU_MU", 0.1, 0.98);
    add_row(rows, 30, "S7_Bsphimumu_0.1_0.98", "DEFAULT", O::S_7_BS_PHI_MU_MU, "S_7_BS_PHI_MU_MU", 0.1, 0.98);
    add_row(rows, 31, "dGamma/dq2_Bsphimumu_1.1_2.5", "DEFAULT", O::DBR_DQ2_BS__PHI_MU_MU, "DBR_DQ2_BS__PHI_MU_MU", 1.1, 2.5);
    add_row(rows, 32, "dGamma/dq2_Bsphimumu_2.5_4", "DEFAULT", O::DBR_DQ2_BS__PHI_MU_MU, "DBR_DQ2_BS__PHI_MU_MU", 2.5, 4);
    add_row(rows, 33, "FL_Bsphimumu_1.1_4", "DEFAULT", O::F_L_BS_PHI_MU_MU, "F_L_BS_PHI_MU_MU", 1.1, 4);
    add_row(rows, 34, "S3_Bsphimumu_1.1_4", "DEFAULT", O::S_3_BS_PHI_MU_MU, "S_3_BS_PHI_MU_MU", 1.1, 4);
    add_row(rows, 35, "S4_Bsphimumu_1.1_4", "DEFAULT", O::S_4_BS_PHI_MU_MU, "S_4_BS_PHI_MU_MU", 1.1, 4);
    add_row(rows, 36, "S7_Bsphimumu_1.1_4", "DEFAULT", O::S_7_BS_PHI_MU_MU, "S_7_BS_PHI_MU_MU", 1.1, 4);
    add_row(rows, 37, "dGamma/dq2_Bsphimumu_4_6", "DEFAULT", O::DBR_DQ2_BS__PHI_MU_MU, "DBR_DQ2_BS__PHI_MU_MU", 4, 6);
    add_row(rows, 38, "FL_Bsphimumu_4_6", "DEFAULT", O::F_L_BS_PHI_MU_MU, "F_L_BS_PHI_MU_MU", 4, 6);
    add_row(rows, 39, "S3_Bsphimumu_4_6", "DEFAULT", O::S_3_BS_PHI_MU_MU, "S_3_BS_PHI_MU_MU", 4, 6);
    add_row(rows, 40, "S4_Bsphimumu_4_6", "DEFAULT", O::S_4_BS_PHI_MU_MU, "S_4_BS_PHI_MU_MU", 4, 6);
    add_row(rows, 41, "S7_Bsphimumu_4_6", "DEFAULT", O::S_7_BS_PHI_MU_MU, "S_7_BS_PHI_MU_MU", 4, 6);
    add_row(rows, 42, "dGamma/dq2_Bsphimumu_15_19", "DEFAULT", O::DBR_DQ2_BS__PHI_MU_MU, "DBR_DQ2_BS__PHI_MU_MU", 15, 19);
    add_row(rows, 43, "FL_Bsphimumu_15_18.9", "DEFAULT", O::F_L_BS_PHI_MU_MU, "F_L_BS_PHI_MU_MU", 15, 18.9);
    add_row(rows, 44, "S3_Bsphimumu_15_18.9", "DEFAULT", O::S_3_BS_PHI_MU_MU, "S_3_BS_PHI_MU_MU", 15, 18.9);
    add_row(rows, 45, "S4_Bsphimumu_15_18.9", "DEFAULT", O::S_4_BS_PHI_MU_MU, "S_4_BS_PHI_MU_MU", 15, 18.9);
    add_row(rows, 46, "S7_Bsphimumu_15_18.9", "DEFAULT", O::S_7_BS_PHI_MU_MU, "S_7_BS_PHI_MU_MU", 15, 18.9);
    add_row(rows, 47, "dGamma/dq2_LambdabLambdamumu_15_20", "DEFAULT", O::DBR_DQ2_LAMBDA_B__LAMBDA_MU_MU, "DBR_DQ2_LAMBDA_B__LAMBDA_MU_MU", 15, 20);
    add_row(rows, 48, "AlFB_LambdabLambdamumu_15_20", "DEFAULT", O::A_FB_L_LAMBDA_B__LAMBDA_MU_MU, "A_FB_L_LAMBDA_B__LAMBDA_MU_MU", 15, 20);
    add_row(rows, 49, "AhFB_LambdabLambdamumu_15_20", "DEFAULT", O::A_FB_H_LAMBDA_B__LAMBDA_MU_MU, "A_FB_H_LAMBDA_B__LAMBDA_MU_MU", 15, 20);
    add_row(rows, 50, "AlhFB_LambdabLambdamumu_15_20", "DEFAULT", O::A_FB_LH_LAMBDA_B__LAMBDA_MU_MU, "A_FB_LH_LAMBDA_B__LAMBDA_MU_MU", 15, 20);
    add_row(rows, 51, "FL_LambdabLambdamumu_15_20", "DEFAULT", O::F_L_LAMBDA_B__LAMBDA_MU_MU, "F_L_LAMBDA_B__LAMBDA_MU_MU", 15, 20);
    add_row(rows, 52, "FL_BKstarmumu_0.1_0.98", "DEFAULT", O::F_L_B__KSTAR_MU_MU, "F_L_B__KSTAR_MU_MU", 0.1, 0.98);
    add_row(rows, 53, "AFB_BKstarmumu_0.1_0.98", "DEFAULT", O::A_FB_B__KSTAR_MU_MU, "A_FB_B__KSTAR_MU_MU", 0.1, 0.98);
    add_row(rows, 54, "S3_BKstarmumu_0.1_0.98", "DEFAULT", O::S_3_B__KSTAR_MU_MU, "S_3_B__KSTAR_MU_MU", 0.1, 0.98);
    add_row(rows, 55, "S4_BKstarmumu_0.1_0.98", "DEFAULT", O::S_4_B__KSTAR_MU_MU, "S_4_B__KSTAR_MU_MU", 0.1, 0.98);
    add_row(rows, 56, "S5_BKstarmumu_0.1_0.98", "DEFAULT", O::S_5_B__KSTAR_MU_MU, "S_5_B__KSTAR_MU_MU", 0.1, 0.98);
    add_row(rows, 57, "S7_BKstarmumu_0.1_0.98", "DEFAULT", O::S_7_B__KSTAR_MU_MU, "S_7_B__KSTAR_MU_MU", 0.1, 0.98);
    add_row(rows, 58, "S8_BKstarmumu_0.1_0.98", "DEFAULT", O::S_8_B__KSTAR_MU_MU, "S_8_B__KSTAR_MU_MU", 0.1, 0.98);
    add_row(rows, 59, "S9_BKstarmumu_0.1_0.98", "DEFAULT", O::S_9_B__KSTAR_MU_MU, "S_9_B__KSTAR_MU_MU", 0.1, 0.98);
    add_row(rows, 60, "FL_BKstarmumu_1.1_2.5", "DEFAULT", O::F_L_B__KSTAR_MU_MU, "F_L_B__KSTAR_MU_MU", 1.1, 2.5);
    add_row(rows, 61, "AFB_BKstarmumu_1.1_2.5", "DEFAULT", O::A_FB_B__KSTAR_MU_MU, "A_FB_B__KSTAR_MU_MU", 1.1, 2.5);
    add_row(rows, 62, "S3_BKstarmumu_1.1_2.5", "DEFAULT", O::S_3_B__KSTAR_MU_MU, "S_3_B__KSTAR_MU_MU", 1.1, 2.5);
    add_row(rows, 63, "S4_BKstarmumu_1.1_2.5", "DEFAULT", O::S_4_B__KSTAR_MU_MU, "S_4_B__KSTAR_MU_MU", 1.1, 2.5);
    add_row(rows, 64, "S5_BKstarmumu_1.1_2.5", "DEFAULT", O::S_5_B__KSTAR_MU_MU, "S_5_B__KSTAR_MU_MU", 1.1, 2.5);
    add_row(rows, 65, "S7_BKstarmumu_1.1_2.5", "DEFAULT", O::S_7_B__KSTAR_MU_MU, "S_7_B__KSTAR_MU_MU", 1.1, 2.5);
    add_row(rows, 66, "S8_BKstarmumu_1.1_2.5", "DEFAULT", O::S_8_B__KSTAR_MU_MU, "S_8_B__KSTAR_MU_MU", 1.1, 2.5);
    add_row(rows, 67, "S9_BKstarmumu_1.1_2.5", "DEFAULT", O::S_9_B__KSTAR_MU_MU, "S_9_B__KSTAR_MU_MU", 1.1, 2.5);
    add_row(rows, 68, "FL_BKstarmumu_2.5_4", "DEFAULT", O::F_L_B__KSTAR_MU_MU, "F_L_B__KSTAR_MU_MU", 2.5, 4);
    add_row(rows, 69, "AFB_BKstarmumu_2.5_4", "DEFAULT", O::A_FB_B__KSTAR_MU_MU, "A_FB_B__KSTAR_MU_MU", 2.5, 4);
    add_row(rows, 70, "S3_BKstarmumu_2.5_4", "DEFAULT", O::S_3_B__KSTAR_MU_MU, "S_3_B__KSTAR_MU_MU", 2.5, 4);
    add_row(rows, 71, "S4_BKstarmumu_2.5_4", "DEFAULT", O::S_4_B__KSTAR_MU_MU, "S_4_B__KSTAR_MU_MU", 2.5, 4);
    add_row(rows, 72, "S5_BKstarmumu_2.5_4", "DEFAULT", O::S_5_B__KSTAR_MU_MU, "S_5_B__KSTAR_MU_MU", 2.5, 4);
    add_row(rows, 73, "S7_BKstarmumu_2.5_4", "DEFAULT", O::S_7_B__KSTAR_MU_MU, "S_7_B__KSTAR_MU_MU", 2.5, 4);
    add_row(rows, 74, "S8_BKstarmumu_2.5_4", "DEFAULT", O::S_8_B__KSTAR_MU_MU, "S_8_B__KSTAR_MU_MU", 2.5, 4);
    add_row(rows, 75, "S9_BKstarmumu_2.5_4", "DEFAULT", O::S_9_B__KSTAR_MU_MU, "S_9_B__KSTAR_MU_MU", 2.5, 4);
    add_row(rows, 76, "FL_BKstarmumu_4_6", "DEFAULT", O::F_L_B__KSTAR_MU_MU, "F_L_B__KSTAR_MU_MU", 4, 6);
    add_row(rows, 77, "AFB_BKstarmumu_4_6", "DEFAULT", O::A_FB_B__KSTAR_MU_MU, "A_FB_B__KSTAR_MU_MU", 4, 6);
    add_row(rows, 78, "S3_BKstarmumu_4_6", "DEFAULT", O::S_3_B__KSTAR_MU_MU, "S_3_B__KSTAR_MU_MU", 4, 6);
    add_row(rows, 79, "S4_BKstarmumu_4_6", "DEFAULT", O::S_4_B__KSTAR_MU_MU, "S_4_B__KSTAR_MU_MU", 4, 6);
    add_row(rows, 80, "S5_BKstarmumu_4_6", "DEFAULT", O::S_5_B__KSTAR_MU_MU, "S_5_B__KSTAR_MU_MU", 4, 6);
    add_row(rows, 81, "S7_BKstarmumu_4_6", "DEFAULT", O::S_7_B__KSTAR_MU_MU, "S_7_B__KSTAR_MU_MU", 4, 6);
    add_row(rows, 82, "S8_BKstarmumu_4_6", "DEFAULT", O::S_8_B__KSTAR_MU_MU, "S_8_B__KSTAR_MU_MU", 4, 6);
    add_row(rows, 83, "S9_BKstarmumu_4_6", "DEFAULT", O::S_9_B__KSTAR_MU_MU, "S_9_B__KSTAR_MU_MU", 4, 6);
    add_row(rows, 84, "FL_BKstarmumu_15_17", "DEFAULT", O::F_L_B__KSTAR_MU_MU, "F_L_B__KSTAR_MU_MU", 15, 17);
    add_row(rows, 85, "AFB_BKstarmumu_15_17", "DEFAULT", O::A_FB_B__KSTAR_MU_MU, "A_FB_B__KSTAR_MU_MU", 15, 17);
    add_row(rows, 86, "S3_BKstarmumu_15_17", "DEFAULT", O::S_3_B__KSTAR_MU_MU, "S_3_B__KSTAR_MU_MU", 15, 17);
    add_row(rows, 87, "S4_BKstarmumu_15_17", "DEFAULT", O::S_4_B__KSTAR_MU_MU, "S_4_B__KSTAR_MU_MU", 15, 17);
    add_row(rows, 88, "S5_BKstarmumu_15_17", "DEFAULT", O::S_5_B__KSTAR_MU_MU, "S_5_B__KSTAR_MU_MU", 15, 17);
    add_row(rows, 89, "S7_BKstarmumu_15_17", "DEFAULT", O::S_7_B__KSTAR_MU_MU, "S_7_B__KSTAR_MU_MU", 15, 17);
    add_row(rows, 90, "S8_BKstarmumu_15_17", "DEFAULT", O::S_8_B__KSTAR_MU_MU, "S_8_B__KSTAR_MU_MU", 15, 17);
    add_row(rows, 91, "S9_BKstarmumu_15_17", "DEFAULT", O::S_9_B__KSTAR_MU_MU, "S_9_B__KSTAR_MU_MU", 15, 17);
    add_row(rows, 92, "FL_BKstarmumu_17_19", "DEFAULT", O::F_L_B__KSTAR_MU_MU, "F_L_B__KSTAR_MU_MU", 17, 19);
    add_row(rows, 93, "AFB_BKstarmumu_17_19", "DEFAULT", O::A_FB_B__KSTAR_MU_MU, "A_FB_B__KSTAR_MU_MU", 17, 19);
    add_row(rows, 94, "S3_BKstarmumu_17_19", "DEFAULT", O::S_3_B__KSTAR_MU_MU, "S_3_B__KSTAR_MU_MU", 17, 19);
    add_row(rows, 95, "S4_BKstarmumu_17_19", "DEFAULT", O::S_4_B__KSTAR_MU_MU, "S_4_B__KSTAR_MU_MU", 17, 19);
    add_row(rows, 96, "S5_BKstarmumu_17_19", "DEFAULT", O::S_5_B__KSTAR_MU_MU, "S_5_B__KSTAR_MU_MU", 17, 19);
    add_row(rows, 97, "S7_BKstarmumu_17_19", "DEFAULT", O::S_7_B__KSTAR_MU_MU, "S_7_B__KSTAR_MU_MU", 17, 19);
    add_row(rows, 98, "S8_BKstarmumu_17_19", "DEFAULT", O::S_8_B__KSTAR_MU_MU, "S_8_B__KSTAR_MU_MU", 17, 19);
    add_row(rows, 99, "S9_BKstarmumu_17_19", "DEFAULT", O::S_9_B__KSTAR_MU_MU, "S_9_B__KSTAR_MU_MU", 17, 19);
    add_row(rows, 100, "R-1_BKstarll_0.045_6", "DEFAULT", O::R_1_B__KSTAR_L_L, "R_1_B__KSTAR_L_L", 0.045, 6);
    add_row(rows, 101, "R-1_B0K0ll_1.1_6", "DEFAULT", O::R_1_B0__K0_L_L, "R_1_B0__K0_L_L", 1.1, 6);
    add_row(rows, 102, "R-1_BKll_1_6_Belle", "Belle", O::R_1_B__K_L_L, "R_1_B__K_L_L", 1, 6);
    add_row(rows, 103, "FH_BKmumu_1_6_CMS", "CMS", O::F_H_B__K_MU_MU, "F_H_B__K_MU_MU", 1, 6);
    add_row(rows, 104, "dGamma/dq2_B0Kstar0ee_0.0009_1", "DEFAULT", O::DBR_DQ2_B0__KSTAR0_E_E, "DBR_DQ2_B0__KSTAR0_E_E", 0.0009, 1);
    add_row(rows, 105, "FL_B0Kstar0ee_0.0008_0.257", "DEFAULT", O::F_L_B0__KSTAR0_E_E, "F_L_B0__KSTAR0_E_E", 0.0008, 0.257);
    add_row(rows, 106, "ATRe_B0Kstar0ee_0.0008_0.257", "DEFAULT", O::A_T_RE_B0__KSTAR0_E_E, "A_T_RE_B0__KSTAR0_E_E", 0.0008, 0.257);
    add_row(rows, 107, "AT2_B0Kstar0ee_0.0008_0.257", "DEFAULT", O::A_T_2_B0__KSTAR0_E_E, "A_T_2_B0__KSTAR0_E_E", 0.0008, 0.257);
    add_row(rows, 108, "dGamma/dq2_BKmumu_0.1_0.98_CMS", "CMS", O::DBR_DQ2_B__K_MU_MU, "DBR_DQ2_B__K_MU_MU", 0.1, 0.98);
    add_row(rows, 109, "dGamma/dq2_BKmumu_1.1_2_CMS", "CMS", O::DBR_DQ2_B__K_MU_MU, "DBR_DQ2_B__K_MU_MU", 1.1, 2);
    add_row(rows, 110, "dGamma/dq2_BKmumu_2_3_CMS", "CMS", O::DBR_DQ2_B__K_MU_MU, "DBR_DQ2_B__K_MU_MU", 2, 3);
    add_row(rows, 111, "dGamma/dq2_BKmumu_3_4_CMS", "CMS", O::DBR_DQ2_B__K_MU_MU, "DBR_DQ2_B__K_MU_MU", 3, 4);
    add_row(rows, 112, "dGamma/dq2_BKmumu_4_5_CMS", "CMS", O::DBR_DQ2_B__K_MU_MU, "DBR_DQ2_B__K_MU_MU", 4, 5);
    add_row(rows, 113, "dGamma/dq2_BKmumu_5_6_CMS", "CMS", O::DBR_DQ2_B__K_MU_MU, "DBR_DQ2_B__K_MU_MU", 5, 6);
    add_row(rows, 114, "dGamma/dq2_BKmumu_14.82_16_CMS", "CMS", O::DBR_DQ2_B__K_MU_MU, "DBR_DQ2_B__K_MU_MU", 14.82, 16);
    add_row(rows, 115, "dGamma/dq2_BKmumu_16_17_CMS", "CMS", O::DBR_DQ2_B__K_MU_MU, "DBR_DQ2_B__K_MU_MU", 16, 17);
    add_row(rows, 116, "dGamma/dq2_BKmumu_17_18_CMS", "CMS", O::DBR_DQ2_B__K_MU_MU, "DBR_DQ2_B__K_MU_MU", 17, 18);
    add_row(rows, 117, "dGamma/dq2_BKmumu_18_19.24_CMS", "CMS", O::DBR_DQ2_B__K_MU_MU, "DBR_DQ2_B__K_MU_MU", 18, 19.24);
    add_row(rows, 118, "dGamma/dq2_BKmumu_19.24_22.9_CMS", "CMS", O::DBR_DQ2_B__K_MU_MU, "DBR_DQ2_B__K_MU_MU", 19.24, 22.9);
    add_row(rows, 119, "R-1_BKll_1.1_6_CMS", "CMS", O::R_1_B__K_L_L, "R_1_B__K_L_L", 1.1, 6);
    add_row(rows, 120, "FL_B0Kstar0mumu_1.1_2_CMS", "CMS", O::F_L_B0__KSTAR0_MU_MU, "F_L_B0__KSTAR0_MU_MU", 1.1, 2);
    add_row(rows, 121, "P1_B0Kstar0mumu_1.1_2_CMS", "CMS", O::P_1_B0__KSTAR0_MU_MU, "P_1_B0__KSTAR0_MU_MU", 1.1, 2);
    add_row(rows, 122, "P2_B0Kstar0mumu_1.1_2_CMS", "CMS", O::P_2_B0__KSTAR0_MU_MU, "P_2_B0__KSTAR0_MU_MU", 1.1, 2);
    add_row(rows, 123, "P3_B0Kstar0mumu_1.1_2_CMS", "CMS", O::P_3_B0__KSTAR0_MU_MU, "P_3_B0__KSTAR0_MU_MU", 1.1, 2);
    add_row(rows, 124, "P4prime_B0Kstar0mumu_1.1_2_CMS", "CMS", O::P_PRIME_4_B0__KSTAR0_MU_MU, "P_PRIME_4_B0__KSTAR0_MU_MU", 1.1, 2);
    add_row(rows, 125, "P5prime_B0Kstar0mumu_1.1_2_CMS", "CMS", O::P_PRIME_5_B0__KSTAR0_MU_MU, "P_PRIME_5_B0__KSTAR0_MU_MU", 1.1, 2);
    add_row(rows, 126, "P6prime_B0Kstar0mumu_1.1_2_CMS", "CMS", O::P_PRIME_6_B0__KSTAR0_MU_MU, "P_PRIME_6_B0__KSTAR0_MU_MU", 1.1, 2);
    add_row(rows, 127, "P8prime_B0Kstar0mumu_1.1_2_CMS", "CMS", O::P_PRIME_8_B0__KSTAR0_MU_MU, "P_PRIME_8_B0__KSTAR0_MU_MU", 1.1, 2);
    add_row(rows, 128, "FL_B0Kstar0mumu_2_4.3_CMS", "CMS", O::F_L_B0__KSTAR0_MU_MU, "F_L_B0__KSTAR0_MU_MU", 2, 4.3);
    add_row(rows, 129, "P1_B0Kstar0mumu_2_4.3_CMS", "CMS", O::P_1_B0__KSTAR0_MU_MU, "P_1_B0__KSTAR0_MU_MU", 2, 4.3);
    add_row(rows, 130, "P2_B0Kstar0mumu_2_4.3_CMS", "CMS", O::P_2_B0__KSTAR0_MU_MU, "P_2_B0__KSTAR0_MU_MU", 2, 4.3);
    add_row(rows, 131, "P3_B0Kstar0mumu_2_4.3_CMS", "CMS", O::P_3_B0__KSTAR0_MU_MU, "P_3_B0__KSTAR0_MU_MU", 2, 4.3);
    add_row(rows, 132, "P4prime_B0Kstar0mumu_2_4.3_CMS", "CMS", O::P_PRIME_4_B0__KSTAR0_MU_MU, "P_PRIME_4_B0__KSTAR0_MU_MU", 2, 4.3);
    add_row(rows, 133, "P5prime_B0Kstar0mumu_2_4.3_CMS", "CMS", O::P_PRIME_5_B0__KSTAR0_MU_MU, "P_PRIME_5_B0__KSTAR0_MU_MU", 2, 4.3);
    add_row(rows, 134, "P6prime_B0Kstar0mumu_2_4.3_CMS", "CMS", O::P_PRIME_6_B0__KSTAR0_MU_MU, "P_PRIME_6_B0__KSTAR0_MU_MU", 2, 4.3);
    add_row(rows, 135, "P8prime_B0Kstar0mumu_2_4.3_CMS", "CMS", O::P_PRIME_8_B0__KSTAR0_MU_MU, "P_PRIME_8_B0__KSTAR0_MU_MU", 2, 4.3);
    add_row(rows, 136, "FL_B0Kstar0mumu_4.3_6_CMS", "CMS", O::F_L_B0__KSTAR0_MU_MU, "F_L_B0__KSTAR0_MU_MU", 4.3, 6);
    add_row(rows, 137, "P1_B0Kstar0mumu_4.3_6_CMS", "CMS", O::P_1_B0__KSTAR0_MU_MU, "P_1_B0__KSTAR0_MU_MU", 4.3, 6);
    add_row(rows, 138, "P2_B0Kstar0mumu_4.3_6_CMS", "CMS", O::P_2_B0__KSTAR0_MU_MU, "P_2_B0__KSTAR0_MU_MU", 4.3, 6);
    add_row(rows, 139, "P3_B0Kstar0mumu_4.3_6_CMS", "CMS", O::P_3_B0__KSTAR0_MU_MU, "P_3_B0__KSTAR0_MU_MU", 4.3, 6);
    add_row(rows, 140, "P4prime_B0Kstar0mumu_4.3_6_CMS", "CMS", O::P_PRIME_4_B0__KSTAR0_MU_MU, "P_PRIME_4_B0__KSTAR0_MU_MU", 4.3, 6);
    add_row(rows, 141, "P5prime_B0Kstar0mumu_4.3_6_CMS", "CMS", O::P_PRIME_5_B0__KSTAR0_MU_MU, "P_PRIME_5_B0__KSTAR0_MU_MU", 4.3, 6);
    add_row(rows, 142, "P6prime_B0Kstar0mumu_4.3_6_CMS", "CMS", O::P_PRIME_6_B0__KSTAR0_MU_MU, "P_PRIME_6_B0__KSTAR0_MU_MU", 4.3, 6);
    add_row(rows, 143, "P8prime_B0Kstar0mumu_4.3_6_CMS", "CMS", O::P_PRIME_8_B0__KSTAR0_MU_MU, "P_PRIME_8_B0__KSTAR0_MU_MU", 4.3, 6);
    add_row(rows, 144, "FL_B0Kstar0mumu_14.18_16_CMS", "CMS", O::F_L_B0__KSTAR0_MU_MU, "F_L_B0__KSTAR0_MU_MU", 14.18, 16);
    add_row(rows, 145, "P1_B0Kstar0mumu_14.18_16_CMS", "CMS", O::P_1_B0__KSTAR0_MU_MU, "P_1_B0__KSTAR0_MU_MU", 14.18, 16);
    add_row(rows, 146, "P2_B0Kstar0mumu_14.18_16_CMS", "CMS", O::P_2_B0__KSTAR0_MU_MU, "P_2_B0__KSTAR0_MU_MU", 14.18, 16);
    add_row(rows, 147, "P3_B0Kstar0mumu_14.18_16_CMS", "CMS", O::P_3_B0__KSTAR0_MU_MU, "P_3_B0__KSTAR0_MU_MU", 14.18, 16);
    add_row(rows, 148, "P4prime_B0Kstar0mumu_14.18_16_CMS", "CMS", O::P_PRIME_4_B0__KSTAR0_MU_MU, "P_PRIME_4_B0__KSTAR0_MU_MU", 14.18, 16);
    add_row(rows, 149, "P5prime_B0Kstar0mumu_14.18_16_CMS", "CMS", O::P_PRIME_5_B0__KSTAR0_MU_MU, "P_PRIME_5_B0__KSTAR0_MU_MU", 14.18, 16);
    add_row(rows, 150, "P6prime_B0Kstar0mumu_14.18_16_CMS", "CMS", O::P_PRIME_6_B0__KSTAR0_MU_MU, "P_PRIME_6_B0__KSTAR0_MU_MU", 14.18, 16);
    add_row(rows, 151, "P8prime_B0Kstar0mumu_14.18_16_CMS", "CMS", O::P_PRIME_8_B0__KSTAR0_MU_MU, "P_PRIME_8_B0__KSTAR0_MU_MU", 14.18, 16);
    add_row(rows, 152, "dGamma/dq2_Bsphiee_0.1_1.1", "DEFAULT", O::DBR_DQ2_BS__PHI_E_E, "DBR_DQ2_BS__PHI_E_E", 0.1, 1.1);
    add_row(rows, 153, "dGamma/dq2_Bsphiee_1.1_6", "DEFAULT", O::DBR_DQ2_BS__PHI_E_E, "DBR_DQ2_BS__PHI_E_E", 1.1, 6);
    add_row(rows, 154, "dGamma/dq2_Bsphiee_15_19", "DEFAULT", O::DBR_DQ2_BS__PHI_E_E, "DBR_DQ2_BS__PHI_E_E", 15, 19);
    add_row(rows, 155, "R-1_Bsphill_0.1_1.1", "DEFAULT", O::R_1_BS__PHI_L_L, "R_1_BS__PHI_L_L", 0.1, 1.1);
    add_row(rows, 156, "R-1_Bsphill_1.1_6", "DEFAULT", O::R_1_BS__PHI_L_L, "R_1_BS__PHI_L_L", 1.1, 6);
    add_row(rows, 157, "R-1_Bsphill_15_19", "DEFAULT", O::R_1_BS__PHI_L_L, "R_1_BS__PHI_L_L", 15, 19);
    add_row(rows, 158, "FL_Bsphiee_0.0009_0.2615", "DEFAULT", O::F_L_BS_PHI_E_E, "F_L_BS_PHI_E_E", 0.0009, 0.2615);
    add_row(rows, 159, "AT2_Bsphiee_0.0009_0.2615", "DEFAULT", O::A_T_2_BS_PHI_E_E, "A_T_2_BS_PHI_E_E", 0.0009, 0.2615);
    add_row(rows, 160, "AT2_B0Kstar0ee_0.0008_1.12_Belle", "Belle", O::A_T_2_B0__KSTAR0_E_E, "A_T_2_B0__KSTAR0_E_E", 0.0008, 1.12);
    add_row(rows, 161, "FL_B0Kstar0ee_1.1_6", "DEFAULT", O::F_L_B0__KSTAR0_E_E, "F_L_B0__KSTAR0_E_E", 1.1, 6);
    add_row(rows, 162, "P1_B0Kstar0ee_1.1_6", "DEFAULT", O::P_1_B0__KSTAR0_E_E, "P_1_B0__KSTAR0_E_E", 1.1, 6);
    add_row(rows, 163, "P2_B0Kstar0ee_1.1_6", "DEFAULT", O::P_2_B0__KSTAR0_E_E, "P_2_B0__KSTAR0_E_E", 1.1, 6);
    add_row(rows, 164, "P3_B0Kstar0ee_1.1_6", "DEFAULT", O::P_3_B0__KSTAR0_E_E, "P_3_B0__KSTAR0_E_E", 1.1, 6);
    add_row(rows, 165, "P4prime_B0Kstar0ee_1.1_6", "DEFAULT", O::P_PRIME_4_B0__KSTAR0_E_E, "P_PRIME_4_B0__KSTAR0_E_E", 1.1, 6);
    add_row(rows, 166, "P5prime_B0Kstar0ee_1.1_6", "DEFAULT", O::P_PRIME_5_B0__KSTAR0_E_E, "P_PRIME_5_B0__KSTAR0_E_E", 1.1, 6);
    add_row(rows, 167, "P6prime_B0Kstar0ee_1.1_6", "DEFAULT", O::P_PRIME_6_B0__KSTAR0_E_E, "P_PRIME_6_B0__KSTAR0_E_E", 1.1, 6);
    add_row(rows, 168, "P8prime_B0Kstar0ee_1.1_6", "DEFAULT", O::P_PRIME_8_B0__KSTAR0_E_E, "P_PRIME_8_B0__KSTAR0_E_E", 1.1, 6);
    add_row(rows, 169, "FL_B0Kstar0mumu_0.06_0.98_LHCb2025c2", "LHCb2025c2", O::F_L_B0__KSTAR0_MU_MU, "F_L_B0__KSTAR0_MU_MU", 0.06, 0.98);
    add_row(rows, 170, "S2s_B0Kstar0mumu_0.06_0.98_LHCb2025c2", "LHCb2025c2", O::S_2S_B0__KSTAR0_MU_MU, "S_2S_B0__KSTAR0_MU_MU", 0.06, 0.98);
    add_row(rows, 171, "S1c_B0Kstar0mumu_0.06_0.98_LHCb2025c2", "LHCb2025c2", O::S_1C_B0__KSTAR0_MU_MU, "S_1C_B0__KSTAR0_MU_MU", 0.06, 0.98);
    add_row(rows, 172, "P1_B0Kstar0mumu_0.06_0.98_LHCb2025c2", "LHCb2025c2", O::P_1_B0__KSTAR0_MU_MU, "P_1_B0__KSTAR0_MU_MU", 0.06, 0.98);
    add_row(rows, 173, "P2_B0Kstar0mumu_0.06_0.98_LHCb2025c2", "LHCb2025c2", O::P_2_B0__KSTAR0_MU_MU, "P_2_B0__KSTAR0_MU_MU", 0.06, 0.98);
    add_row(rows, 174, "P3_B0Kstar0mumu_0.06_0.98_LHCb2025c2", "LHCb2025c2", O::P_3_B0__KSTAR0_MU_MU, "P_3_B0__KSTAR0_MU_MU", 0.06, 0.98);
    add_row(rows, 175, "P4prime_B0Kstar0mumu_0.06_0.98_LHCb2025c2", "LHCb2025c2", O::P_PRIME_4_B0__KSTAR0_MU_MU, "P_PRIME_4_B0__KSTAR0_MU_MU", 0.06, 0.98);
    add_row(rows, 176, "P5prime_B0Kstar0mumu_0.06_0.98_LHCb2025c2", "LHCb2025c2", O::P_PRIME_5_B0__KSTAR0_MU_MU, "P_PRIME_5_B0__KSTAR0_MU_MU", 0.06, 0.98);
    add_row(rows, 177, "P6prime_B0Kstar0mumu_0.06_0.98_LHCb2025c2", "LHCb2025c2", O::P_PRIME_6_B0__KSTAR0_MU_MU, "P_PRIME_6_B0__KSTAR0_MU_MU", 0.06, 0.98);
    add_row(rows, 178, "P8prime_B0Kstar0mumu_0.06_0.98_LHCb2025c2", "LHCb2025c2", O::P_PRIME_8_B0__KSTAR0_MU_MU, "P_PRIME_8_B0__KSTAR0_MU_MU", 0.06, 0.98);
    add_row(rows, 179, "S6c_B0Kstar0mumu_0.06_0.98_LHCb2025c2", "LHCb2025c2", O::S_6C_B0__KSTAR0_MU_MU, "S_6C_B0__KSTAR0_MU_MU", 0.06, 0.98);
    add_row(rows, 180, "dGamma/dq2_B0Kstar0mumu_0.06_0.98_LHCb2025c2", "LHCb2025c2", O::DBR_DQ2_B0__KSTAR0_MU_MU, "DBR_DQ2_B0__KSTAR0_MU_MU", 0.06, 0.98);
    add_row(rows, 181, "FL_B0Kstar0mumu_1.1_2.5_LHCb2025c2", "LHCb2025c2", O::F_L_B0__KSTAR0_MU_MU, "F_L_B0__KSTAR0_MU_MU", 1.1, 2.5);
    add_row(rows, 182, "S1c_B0Kstar0mumu_1.1_2.5_LHCb2025c2", "LHCb2025c2", O::S_1C_B0__KSTAR0_MU_MU, "S_1C_B0__KSTAR0_MU_MU", 1.1, 2.5);
    add_row(rows, 183, "P1_B0Kstar0mumu_1.1_2.5_LHCb2025c2", "LHCb2025c2", O::P_1_B0__KSTAR0_MU_MU, "P_1_B0__KSTAR0_MU_MU", 1.1, 2.5);
    add_row(rows, 184, "P2_B0Kstar0mumu_1.1_2.5_LHCb2025c2", "LHCb2025c2", O::P_2_B0__KSTAR0_MU_MU, "P_2_B0__KSTAR0_MU_MU", 1.1, 2.5);
    add_row(rows, 185, "P3_B0Kstar0mumu_1.1_2.5_LHCb2025c2", "LHCb2025c2", O::P_3_B0__KSTAR0_MU_MU, "P_3_B0__KSTAR0_MU_MU", 1.1, 2.5);
    add_row(rows, 186, "P4prime_B0Kstar0mumu_1.1_2.5_LHCb2025c2", "LHCb2025c2", O::P_PRIME_4_B0__KSTAR0_MU_MU, "P_PRIME_4_B0__KSTAR0_MU_MU", 1.1, 2.5);
    add_row(rows, 187, "P5prime_B0Kstar0mumu_1.1_2.5_LHCb2025c2", "LHCb2025c2", O::P_PRIME_5_B0__KSTAR0_MU_MU, "P_PRIME_5_B0__KSTAR0_MU_MU", 1.1, 2.5);
    add_row(rows, 188, "P6prime_B0Kstar0mumu_1.1_2.5_LHCb2025c2", "LHCb2025c2", O::P_PRIME_6_B0__KSTAR0_MU_MU, "P_PRIME_6_B0__KSTAR0_MU_MU", 1.1, 2.5);
    add_row(rows, 189, "P8prime_B0Kstar0mumu_1.1_2.5_LHCb2025c2", "LHCb2025c2", O::P_PRIME_8_B0__KSTAR0_MU_MU, "P_PRIME_8_B0__KSTAR0_MU_MU", 1.1, 2.5);
    add_row(rows, 190, "dGamma/dq2_B0Kstar0mumu_1.1_2.5_LHCb2025c2", "LHCb2025c2", O::DBR_DQ2_B0__KSTAR0_MU_MU, "DBR_DQ2_B0__KSTAR0_MU_MU", 1.1, 2.5);
    add_row(rows, 191, "FL_B0Kstar0mumu_2.5_4.0_LHCb2025c2", "LHCb2025c2", O::F_L_B0__KSTAR0_MU_MU, "F_L_B0__KSTAR0_MU_MU", 2.5, 4);
    add_row(rows, 192, "S1c_B0Kstar0mumu_2.5_4.0_LHCb2025c2", "LHCb2025c2", O::S_1C_B0__KSTAR0_MU_MU, "S_1C_B0__KSTAR0_MU_MU", 2.5, 4);
    add_row(rows, 193, "P1_B0Kstar0mumu_2.5_4.0_LHCb2025c2", "LHCb2025c2", O::P_1_B0__KSTAR0_MU_MU, "P_1_B0__KSTAR0_MU_MU", 2.5, 4);
    add_row(rows, 194, "P2_B0Kstar0mumu_2.5_4.0_LHCb2025c2", "LHCb2025c2", O::P_2_B0__KSTAR0_MU_MU, "P_2_B0__KSTAR0_MU_MU", 2.5, 4);
    add_row(rows, 195, "P3_B0Kstar0mumu_2.5_4.0_LHCb2025c2", "LHCb2025c2", O::P_3_B0__KSTAR0_MU_MU, "P_3_B0__KSTAR0_MU_MU", 2.5, 4);
    add_row(rows, 196, "P4prime_B0Kstar0mumu_2.5_4.0_LHCb2025c2", "LHCb2025c2", O::P_PRIME_4_B0__KSTAR0_MU_MU, "P_PRIME_4_B0__KSTAR0_MU_MU", 2.5, 4);
    add_row(rows, 197, "P5prime_B0Kstar0mumu_2.5_4.0_LHCb2025c2", "LHCb2025c2", O::P_PRIME_5_B0__KSTAR0_MU_MU, "P_PRIME_5_B0__KSTAR0_MU_MU", 2.5, 4);
    add_row(rows, 198, "P6prime_B0Kstar0mumu_2.5_4.0_LHCb2025c2", "LHCb2025c2", O::P_PRIME_6_B0__KSTAR0_MU_MU, "P_PRIME_6_B0__KSTAR0_MU_MU", 2.5, 4);
    add_row(rows, 199, "P8prime_B0Kstar0mumu_2.5_4.0_LHCb2025c2", "LHCb2025c2", O::P_PRIME_8_B0__KSTAR0_MU_MU, "P_PRIME_8_B0__KSTAR0_MU_MU", 2.5, 4);
    add_row(rows, 200, "dGamma/dq2_B0Kstar0mumu_2.5_4.0_LHCb2025c2", "LHCb2025c2", O::DBR_DQ2_B0__KSTAR0_MU_MU, "DBR_DQ2_B0__KSTAR0_MU_MU", 2.5, 4);
    add_row(rows, 201, "FL_B0Kstar0mumu_4.0_6.0_LHCb2025c2", "LHCb2025c2", O::F_L_B0__KSTAR0_MU_MU, "F_L_B0__KSTAR0_MU_MU", 4, 6);
    add_row(rows, 202, "S1c_B0Kstar0mumu_4.0_6.0_LHCb2025c2", "LHCb2025c2", O::S_1C_B0__KSTAR0_MU_MU, "S_1C_B0__KSTAR0_MU_MU", 4, 6);
    add_row(rows, 203, "P1_B0Kstar0mumu_4.0_6.0_LHCb2025c2", "LHCb2025c2", O::P_1_B0__KSTAR0_MU_MU, "P_1_B0__KSTAR0_MU_MU", 4, 6);
    add_row(rows, 204, "P2_B0Kstar0mumu_4.0_6.0_LHCb2025c2", "LHCb2025c2", O::P_2_B0__KSTAR0_MU_MU, "P_2_B0__KSTAR0_MU_MU", 4, 6);
    add_row(rows, 205, "P3_B0Kstar0mumu_4.0_6.0_LHCb2025c2", "LHCb2025c2", O::P_3_B0__KSTAR0_MU_MU, "P_3_B0__KSTAR0_MU_MU", 4, 6);
    add_row(rows, 206, "P4prime_B0Kstar0mumu_4.0_6.0_LHCb2025c2", "LHCb2025c2", O::P_PRIME_4_B0__KSTAR0_MU_MU, "P_PRIME_4_B0__KSTAR0_MU_MU", 4, 6);
    add_row(rows, 207, "P5prime_B0Kstar0mumu_4.0_6.0_LHCb2025c2", "LHCb2025c2", O::P_PRIME_5_B0__KSTAR0_MU_MU, "P_PRIME_5_B0__KSTAR0_MU_MU", 4, 6);
    add_row(rows, 208, "P6prime_B0Kstar0mumu_4.0_6.0_LHCb2025c2", "LHCb2025c2", O::P_PRIME_6_B0__KSTAR0_MU_MU, "P_PRIME_6_B0__KSTAR0_MU_MU", 4, 6);
    add_row(rows, 209, "P8prime_B0Kstar0mumu_4.0_6.0_LHCb2025c2", "LHCb2025c2", O::P_PRIME_8_B0__KSTAR0_MU_MU, "P_PRIME_8_B0__KSTAR0_MU_MU", 4, 6);
    add_row(rows, 210, "dGamma/dq2_B0Kstar0mumu_4.0_6.0_LHCb2025c2", "LHCb2025c2", O::DBR_DQ2_B0__KSTAR0_MU_MU, "DBR_DQ2_B0__KSTAR0_MU_MU", 4, 6);
    add_row(rows, 211, "FL_B0Kstar0mumu_15.0_17.0_LHCb2025c2", "LHCb2025c2", O::F_L_B0__KSTAR0_MU_MU, "F_L_B0__KSTAR0_MU_MU", 15, 17);
    add_row(rows, 212, "S1c_B0Kstar0mumu_15.0_17.0_LHCb2025c2", "LHCb2025c2", O::S_1C_B0__KSTAR0_MU_MU, "S_1C_B0__KSTAR0_MU_MU", 15, 17);
    add_row(rows, 213, "P1_B0Kstar0mumu_15.0_17.0_LHCb2025c2", "LHCb2025c2", O::P_1_B0__KSTAR0_MU_MU, "P_1_B0__KSTAR0_MU_MU", 15, 17);
    add_row(rows, 214, "P2_B0Kstar0mumu_15.0_17.0_LHCb2025c2", "LHCb2025c2", O::P_2_B0__KSTAR0_MU_MU, "P_2_B0__KSTAR0_MU_MU", 15, 17);
    add_row(rows, 215, "P3_B0Kstar0mumu_15.0_17.0_LHCb2025c2", "LHCb2025c2", O::P_3_B0__KSTAR0_MU_MU, "P_3_B0__KSTAR0_MU_MU", 15, 17);
    add_row(rows, 216, "P4prime_B0Kstar0mumu_15.0_17.0_LHCb2025c2", "LHCb2025c2", O::P_PRIME_4_B0__KSTAR0_MU_MU, "P_PRIME_4_B0__KSTAR0_MU_MU", 15, 17);
    add_row(rows, 217, "P5prime_B0Kstar0mumu_15.0_17.0_LHCb2025c2", "LHCb2025c2", O::P_PRIME_5_B0__KSTAR0_MU_MU, "P_PRIME_5_B0__KSTAR0_MU_MU", 15, 17);
    add_row(rows, 218, "P6prime_B0Kstar0mumu_15.0_17.0_LHCb2025c2", "LHCb2025c2", O::P_PRIME_6_B0__KSTAR0_MU_MU, "P_PRIME_6_B0__KSTAR0_MU_MU", 15, 17);
    add_row(rows, 219, "P8prime_B0Kstar0mumu_15.0_17.0_LHCb2025c2", "LHCb2025c2", O::P_PRIME_8_B0__KSTAR0_MU_MU, "P_PRIME_8_B0__KSTAR0_MU_MU", 15, 17);
    add_row(rows, 220, "dGamma/dq2_B0Kstar0mumu_15.0_17.0_LHCb2025c2", "LHCb2025c2", O::DBR_DQ2_B0__KSTAR0_MU_MU, "DBR_DQ2_B0__KSTAR0_MU_MU", 15, 17);
    add_row(rows, 221, "FL_B0Kstar0mumu_17.0_19.0_LHCb2025c2", "LHCb2025c2", O::F_L_B0__KSTAR0_MU_MU, "F_L_B0__KSTAR0_MU_MU", 17, 19);
    add_row(rows, 222, "S1c_B0Kstar0mumu_17.0_19.0_LHCb2025c2", "LHCb2025c2", O::S_1C_B0__KSTAR0_MU_MU, "S_1C_B0__KSTAR0_MU_MU", 17, 19);
    add_row(rows, 223, "P1_B0Kstar0mumu_17.0_19.0_LHCb2025c2", "LHCb2025c2", O::P_1_B0__KSTAR0_MU_MU, "P_1_B0__KSTAR0_MU_MU", 17, 19);
    add_row(rows, 224, "P2_B0Kstar0mumu_17.0_19.0_LHCb2025c2", "LHCb2025c2", O::P_2_B0__KSTAR0_MU_MU, "P_2_B0__KSTAR0_MU_MU", 17, 19);
    add_row(rows, 225, "P3_B0Kstar0mumu_17.0_19.0_LHCb2025c2", "LHCb2025c2", O::P_3_B0__KSTAR0_MU_MU, "P_3_B0__KSTAR0_MU_MU", 17, 19);
    add_row(rows, 226, "P4prime_B0Kstar0mumu_17.0_19.0_LHCb2025c2", "LHCb2025c2", O::P_PRIME_4_B0__KSTAR0_MU_MU, "P_PRIME_4_B0__KSTAR0_MU_MU", 17, 19);
    add_row(rows, 227, "P5prime_B0Kstar0mumu_17.0_19.0_LHCb2025c2", "LHCb2025c2", O::P_PRIME_5_B0__KSTAR0_MU_MU, "P_PRIME_5_B0__KSTAR0_MU_MU", 17, 19);
    add_row(rows, 228, "P6prime_B0Kstar0mumu_17.0_19.0_LHCb2025c2", "LHCb2025c2", O::P_PRIME_6_B0__KSTAR0_MU_MU, "P_PRIME_6_B0__KSTAR0_MU_MU", 17, 19);
    add_row(rows, 229, "P8prime_B0Kstar0mumu_17.0_19.0_LHCb2025c2", "LHCb2025c2", O::P_PRIME_8_B0__KSTAR0_MU_MU, "P_PRIME_8_B0__KSTAR0_MU_MU", 17, 19);
    add_row(rows, 230, "dGamma/dq2_B0Kstar0mumu_17.0_19.0_LHCb2025c2", "LHCb2025c2", O::DBR_DQ2_B0__KSTAR0_MU_MU, "DBR_DQ2_B0__KSTAR0_MU_MU", 17, 19);

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

    const std::string out_path = (argc >= 2) ? argv[1] : "hyperiso_allobs.csv";
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
                        oi.add_observable(
                            BinnedObservableId{ObservableMapper::to_id(r.hyp_obs), {r.q2min, r.q2max}},
                            QCDOrder::NNLO,
                            add_deps
                        );
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

        out << "code,idx,si_name,experiment,hyper_enum,bin_low,bin_high,value,err_down,err_up,dC7,dC8,dC9,dC10,lha,status\n";
        out << std::setprecision(17);

        int nok = 0;
        int nfail = 0;
        for (const auto& r : rows) {
            double value = std::numeric_limits<double>::quiet_NaN();
            std::string status = "ok";
            try {
                value = find_value(all, r);
                ++nok;
            } catch (const std::exception& e) {
                status = e.what();
                ++nfail;
            }

            out << "hyperiso," << r.idx << ","
                << csv_escape(r.si_name) << ","
                << csv_escape(r.experiment) << ","
                << csv_escape(r.hyp_enum) << ",";
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

        std::cout << "CSV ecrit dans " << out_path << " (ok=" << nok << ", fail=" << nfail << ")\n";
        return nfail == 0 ? 0 : 2;
    } catch (const std::exception& e) {
        std::cerr << "[FATAL] " << e.what() << "\n";
        return 1;
    }
}
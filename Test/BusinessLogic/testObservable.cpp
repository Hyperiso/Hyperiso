#include "MemoryManager.h"
#include "Parameters.h"
#include "Logger.h"
#include "ModelEvaluator.h"
#include "BKstarDecay.h"
#include "BllDecay.h"
#include "Observable.h"

#include <iostream>
#include <cassert>

int main() {

    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    auto mm = MemoryManager::GetInstance();  // Initialize program manager with LHA file containing SMINPUTS block
    mm->init("Test/InputFiles/testInput.flha", Model::SM);  // Initialize parameters from given LHA file

    auto sm_p = Parameters::GetInstance(ParameterType::SM);
    auto flavp = Parameters::GetInstance(ParameterType::FLAVOR);
    double m_Bs = (*flavp)("FMASS", 531);
    double m_Bd = (*flavp)("FMASS", 511);
    double m_Bu = (*flavp)("FMASS", 521);

    auto B_Ks = std::make_shared<BKstarDecay>(QCDOrder::NNLO, 80, QCDHelper::mass_b_1S() / 2);
    auto B_ll = std::make_shared<BllDecay>(QCDOrder::NNLO, 80, m_Bd);

    Observable delta_0 (Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, B_Ks);
    Observable br_Bs__mu_mu(Observables::BR_BS_MUMU, B_ll);
    Observable br_Bs__mu_mu__untag(Observables::BR_BS_MUMU_UNTAG, B_ll);
    Observable br_Bd__mu_mu(Observables::BR_BD_MUMU, B_ll);
    // delta_0.add_dependences({
    //     ParamId{ParameterType::FF, "B_Ks", 1},
    //     ParamId{ParameterType::FF, "B_Ks", 2},
    //     ParamId{ParameterType::FF, "B_Ks", 3},
    //     ParamId{ParameterType::FF, "B_Ks", 4},
    //     ParamId{ParameterType::FF, "B_Ks", 10},
    //     ParamId{ParameterType::FF, "B_Ks", 11},
    //     ParamId{ParameterType::FF, "B_Ks", 12},
    //     ParamId{ParameterType::FLAVOR, "FCONST", 32301},
    //     ParamId{ParameterType::FLAVOR, "FCONST", 32302},
    //     ParamId{ParameterType::FLAVOR, "FCONST", 52101},
    //     ParamId{ParameterType::FLAVOR, "FMASS", 323},
    //     ParamId{ParameterType::FLAVOR, "FMASS", 521},
    //     ParamId{ParameterType::SM, "MASS", 4},
    //     ParamId{ParameterType::SM, "MASS", 5},
    // });

    LOG_INFO("---------- Observables -----------");
    // LOG_INFO("BR(Bs > mu mu) = ",           br_Bs__mu_mu.eval(),        " +- ", std::sqrt(br_Bs__mu_mu.variance()));
    // LOG_INFO("Detailed uncertainties:");
    // for (auto &&p : br_Bs__mu_mu.get_leading_uncertainties(3)) {
    //     LOG_INFO("\t(", p.first.first, ",", p.first.second, ") :", p.second);
    // }

    LOG_INFO("Untagged BR(Bs > mu mu) = ",  br_Bs__mu_mu.eval(), " +- ", std::sqrt(br_Bs__mu_mu.variance()));
    LOG_INFO("Untagged BR(Bs > mu mu) = ",  br_Bs__mu_mu__untag.eval(), " +- ", std::sqrt(br_Bs__mu_mu__untag.variance()));
    LOG_INFO("BR(Bd > mu mu) = ",           br_Bd__mu_mu.eval(),        " +- ", std::sqrt(br_Bd__mu_mu.variance()));
    // LOG_INFO("BR(Bu > tau nu) = ",          br_Bu__tau_nu.eval(),       " +- ", std::sqrt(br_Bu__tau_nu.variance()));
    LOG_INFO("Delta_0(B > K* gamma) = ",       delta_0.eval(),             " +- ", std::sqrt(delta_0.variance()));
    // LOG_INFO("Delta_0(B > K* gamma) = ", delta_0.eval());
    // LOG_INFO("Detailed uncertainties:");
    // for (auto &&p : delta_0.get_leading_uncertainties(3)) {
    //     LOG_INFO("\t(", p.first.block, ",", p.first.code, ") :", p.second);
    // }

    // delta_0.print_gradient(std::cout);
}
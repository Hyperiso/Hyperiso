#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>

#include "CorrelationRepo.h"
#include "Include.h"

#include "registry_init.hpp"
#include "observable_ids.hpp"

static ExperimentObs eobs(const std::string& exp, Observables o) {
    return ExperimentObs{exp, o};
}

int main() {
    std::cout << "== Running CorrelationRepository integration tests ==\n";

    init_all_builtins();

    auto repo = std::make_unique<CorrelationRepository>();

    const std::string exp1 = "LHCb";
    const std::string exp2 = "CMS";

    auto eo1      = eobs(exp1, Observables::BR_BS_MUMU);
    auto eo2      = eobs(exp1, Observables::BR_BD_MUMU);
    auto eo3      = eobs(exp1, Observables::R_D);
    auto eo2_exp2 = eobs(exp2, Observables::BR_BD_MUMU);

    auto obsMat = std::make_shared<CorrelationMatrixPair<ExperimentObs>>();
    obsMat->emplace({eo1, eo2}, 0.30, 0.40);
    repo->set_correlation_matrix(obsMat);

    {
        auto v = repo->get_correlation(exp1,
                                       Observables::BR_BS_MUMU,
                                       Observables::BR_BD_MUMU);
        assert(std::abs(v.first  - 0.30) < 1e-12);
        assert(std::abs(v.second - 0.40) < 1e-12);

        double comb = repo->get_combined_correlation(eo1, eo2);
        assert(std::abs(comb - std::hypot(0.30, 0.40)) < 1e-12);
    }

    auto obsMat2 = std::make_shared<CorrelationMatrixPair<ExperimentObs>>();
    obsMat2->emplace({eo1, eo2},      0.95, 0.05); // overwrite existing
    obsMat2->emplace({eo2, eo3},      0.10, 0.20); // same experiment
    obsMat2->emplace({eo2, eo2_exp2}, 0.60, 0.70); // cross-experiment
    repo->merge_correlation_matrix(obsMat2);

    {
        auto v12 = repo->get_correlation(exp1,
                                         Observables::BR_BS_MUMU,
                                         Observables::BR_BD_MUMU);
        assert(std::abs(v12.first  - 0.95) < 1e-12);
        assert(std::abs(v12.second - 0.05) < 1e-12);

        auto v23 = repo->get_correlation(exp1,
                                         Observables::BR_BD_MUMU,
                                         Observables::R_D);
        auto v32 = repo->get_correlation(exp1,
                                         Observables::R_D,
                                         Observables::BR_BD_MUMU);
        assert(std::abs(v23.first  - 0.10) < 1e-12);
        assert(std::abs(v23.second - 0.20) < 1e-12);
        assert(std::abs(v32.first  - 0.10) < 1e-12);
        assert(std::abs(v32.second - 0.20) < 1e-12);

        auto vx = repo->get_correlation(exp1, eo2.obs, exp2, eo2_exp2.obs);
        auto vy = repo->get_correlation(exp2, eo2_exp2.obs, exp1, eo2.obs);
        assert(std::abs(vx.first  - 0.60) < 1e-12);
        assert(std::abs(vx.second - 0.70) < 1e-12);
        assert(std::abs(vy.first  - 0.60) < 1e-12);
        assert(std::abs(vy.second - 0.70) < 1e-12);
    }

    auto p1 = ParamId{ParameterType::SM, "GAUGE", LhaID(1)};
    auto p2 = ParamId{ParameterType::SM, "GAUGE", LhaID(2)};
    auto p3 = ParamId{ParameterType::SM, "MASS",  LhaID(3)};

    auto parMat = std::make_shared<CorrelationMatrixPair<ParamId>>();
    parMat->emplace({p1, p2}, 0.12, 0.16);
    repo->set_correlation_matrix(parMat);

    {
        auto v = repo->get_correlation(p2, p1);
        assert(std::abs(v.first  - 0.12) < 1e-12);
        assert(std::abs(v.second - 0.16) < 1e-12);

        double comb = repo->get_combined_correlation(p1, p2);
        assert(std::abs(comb - std::hypot(0.12, 0.16)) < 1e-12);
    }

    auto parMat2 = std::make_shared<CorrelationMatrixPair<ParamId>>();
    parMat2->emplace({p1, p2}, 0.50, 0.00);
    parMat2->emplace({p2, p3}, 0.25, 0.30);
    repo->merge_correlation_matrix(parMat2);

    {
        auto v12 = repo->get_correlation(p1, p2);
        auto v23 = repo->get_correlation(p2, p3);
        assert(std::abs(v12.first  - 0.50) < 1e-12);
        assert(std::abs(v12.second - 0.00) < 1e-12);
        assert(std::abs(v23.first  - 0.25) < 1e-12);
        assert(std::abs(v23.second - 0.30) < 1e-12);
    }

    repo->print_content();

    std::cout << "\nCorrelationRepository integration suite passed.\n";
    return 0;
}
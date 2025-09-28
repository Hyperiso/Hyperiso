#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>

#include "CorrelationRepo.h"
#include "General.h"

static ParamId pid(const std::string& block, long code){
    return ParamId{ParameterType::SM, block, LhaID(code)};
}
static ObservableId oid(const std::string& name){
    return ObservableId(name); 
}

int main(){
    std::cout << "== Running CorrelationRepository unit tests ==\n";

    CorrelationMatrixPair<ParamId> cmpP;

    auto pA = pid("GAUGE", 1);
    auto pB = pid("GAUGE", 2);
    auto pC = pid("MASS",  3);

    {
        auto v = cmpP.at({pA, pB});
        assert(v.first  == 0.0);
        assert(v.second == 0.0);
    }

    cmpP.emplace({pA, pB}, /*stat=*/0.4, /*syst=*/0.3);
    {
        auto v1 = cmpP.at({pA, pB});
        auto v2 = cmpP.at({pB, pA});
        assert(std::abs(v1.first  - 0.4) < 1e-12 && std::abs(v1.second - 0.3) < 1e-12);
        assert(std::abs(v2.first  - 0.4) < 1e-12 && std::abs(v2.second - 0.3) < 1e-12);
    }

    cmpP.emplace({pA, pA}, 1.0, 0.0);
    {
        auto vd = cmpP.at({pA, pA});
        assert(std::abs(vd.first - 1.0) < 1e-12 && std::abs(vd.second - 0.0) < 1e-12);
    }

    auto repo = std::make_unique<CorrelationRepository>();
    repo->set_correlation_matrix(std::make_shared<CorrelationMatrixPair<ParamId>>(cmpP));

    {
        auto v = repo->get_correlation(pA, pB);
        assert(std::abs(v.first - 0.4) < 1e-12 && std::abs(v.second - 0.3) < 1e-12);

        double comb = repo->get_combined_correlation(pA, pB);
        assert(std::abs(comb - std::hypot(0.4, 0.3)) < 1e-12);
    }


    auto addP = std::make_shared<CorrelationMatrixPair<ParamId>>();
    addP->emplace({pA, pB}, 0.9, 0.1);
    addP->emplace({pB, pC}, 0.2, 0.7);

    repo->merge_correlation_matrix(addP);

    {
        auto vAB = repo->get_correlation(pA, pB);
        assert(std::abs(vAB.first - 0.9) < 1e-12 && std::abs(vAB.second - 0.1) < 1e-12);
        auto vBA = repo->get_correlation(pB, pA);
        assert(std::abs(vBA.first - 0.9) < 1e-12 && std::abs(vBA.second - 0.1) < 1e-12);

        auto vBC = repo->get_correlation(pB, pC);
        auto vCB = repo->get_correlation(pC, pB);
        assert(std::abs(vBC.first - 0.2) < 1e-12 && std::abs(vBC.second - 0.7) < 1e-12);
        assert(std::abs(vCB.first - 0.2) < 1e-12 && std::abs(vCB.second - 0.7) < 1e-12);
    }


    CorrelationMatrixPair<ObservableId> cmpO;

    auto o1 = oid("OBS_A");
    auto o2 = oid("OBS_B");

    auto v0 = cmpO.at({o1, o2});
    assert(v0.first == 0.0 && v0.second == 0.0);

    cmpO.emplace({o1, o2}, 0.33, 0.44);

    auto repoObs = std::make_unique<CorrelationRepository>();
    repoObs->set_correlation_matrix(std::make_shared<CorrelationMatrixPair<ObservableId>>(cmpO));

    {
        auto v = repoObs->get_correlation(o1, o2);
        assert(std::abs(v.first - 0.33) < 1e-12 && std::abs(v.second - 0.44) < 1e-12);

        double comb = repoObs->get_combined_correlation(o1, o2);
        assert(std::abs(comb - std::hypot(0.33, 0.44)) < 1e-12);
    }

    auto addO = std::make_shared<CorrelationMatrixPair<ObservableId>>();
    auto o3 = oid("OBS_C");
    addO->emplace({o1, o2}, 0.7, 0.0);
    addO->emplace({o2, o3}, 0.1, 0.2);
    repoObs->merge_correlation_matrix(addO);

    {
        auto v12 = repoObs->get_correlation(o1, o2);
        assert(std::abs(v12.first - 0.7) < 1e-12 && std::abs(v12.second - 0.0) < 1e-12);

        auto v23 = repoObs->get_correlation(o2, o3);
        auto v32 = repoObs->get_correlation(o3, o2);
        assert(std::abs(v23.first - 0.1) < 1e-12 && std::abs(v23.second - 0.2) < 1e-12);
        assert(std::abs(v32.first - 0.1) < 1e-12 && std::abs(v32.second - 0.2) < 1e-12);
    }

    std::cout << "\n CorrelationRepository unit suite passed.\n";
    return 0;
}

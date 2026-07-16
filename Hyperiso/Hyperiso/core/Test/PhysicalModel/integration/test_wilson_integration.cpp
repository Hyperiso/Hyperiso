#include <cassert>
#include <iostream>
#include <unordered_map>
#include <memory>
#include <cmath>

#include "Wilson.h"  
#include "Parameter.h"
#include "IParameterProxy.h"
#include "Math.h"
#include "BWilson.h" 

static std::shared_ptr<Parameter> P(double x) {
    return std::make_shared<Parameter>(ParamId("DUMMY", 0), x, 0.0, 0.0);
}

class SpyProxy : public IParameterProxy<std::string, LhaID> {
public:
    mutable std::string last_block;
    mutable LhaID last_id{0};
    double ret = 7.0;

    scalar_t operator()(const std::string& blk, const LhaID& id) const override {
        last_block = blk;
        last_id    = id;
        return ret;
    }
    bool exist(const std::string&, const LhaID&) const override { return true; }
    double get_scale(const std::string&) const override { return 0.0; }
};

int main() {
    std::cout << "== C7 INTEGRATION ==\n";

    C7 c7;

    {
        auto sLO   = c7.get_sources(QCDOrder::LO);
        auto sNLO  = c7.get_sources(QCDOrder::NLO);
        auto sNNLO = c7.get_sources(QCDOrder::NNLO);

        assert(sLO.count({ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2, 1)}) == 1);

        assert(sNLO.count({"WPARAM_MATCH_SM", LhaID(3)}) == 1);
        assert(sNLO.count({"WPARAM_MATCH_SM", LhaID(6)}) == 1);
        assert(sNLO.count({"WPARAM_MATCH_SM", LhaID(2,1)}) == 1);
        assert(sNLO.count({"EW_SCALE",        LhaID(1)}) == 1);

        assert(sNNLO.count({"WPARAM_MATCH_SM", LhaID(7)}) == 1);
        assert(sNNLO.count({"WPARAM_MATCH_SM", LhaID(8)}) == 1);
        assert(sNNLO.count({ParameterType::SM,     "MASS",            LhaID(24)}) == 1);
    }

    {
        std::unordered_map<ParamId, std::shared_ptr<Parameter>> src;
        src[{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2,1)}] = P(4.5);

        auto fLO = c7.get_func(QCDOrder::LO);
        double v1 = fLO(ParamSrc(src));
        src[{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2,1)}] = P(5.0);
        double v2 = fLO(ParamSrc(src));
        assert(v1 != v2);
    }

    {
        std::unordered_map<ParamId, std::shared_ptr<Parameter>> src;
        src[{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(3)}]    = P(0.1);
        src[{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(6)}]    = P(173.0);
        src[{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2,1)}]  = P(4.7);
        src[{ParameterType::WILSON, "EW_SCALE",        LhaID(1)}]    = P(80.0);

        auto fNLO = c7.get_func(QCDOrder::NLO);
        double v1 = fNLO(ParamSrc(src));
        src[{ParameterType::WILSON, "EW_SCALE", LhaID(1)}] = P(100.0);
        double v2 = fNLO(ParamSrc(src));
        assert(v1 != v2);
    }

    {
        std::unordered_map<ParamId, std::shared_ptr<Parameter>> src;
        src[{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(3)}]    = P(0.1);
        src[{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(6)}]    = P(173.0);
        src[{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2,1)}]  = P(4.7);
        src[{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(7)}]    = P((173.0*173.0)/(80.379*80.379));
        src[{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(8)}]    = P(1.0);
        src[{ParameterType::WILSON, "EW_SCALE",        LhaID(1)}]    = P(80.0);
        src[{ParameterType::SM,     "MASS",            LhaID(24)}]   = P(80.379);

        auto fNNLO = c7.get_func(QCDOrder::NNLO);
        double v1 = fNNLO(ParamSrc(src));
        src[{ParameterType::SM, "MASS", LhaID(24)}] = P(82.0);
        double v2 = fNNLO(ParamSrc(src));
        assert(v1 != v2);
    }

    {
        assert(c7.get_lhaid(QCDOrder::LO)   == LhaID(305, 4422, 0, 0));
        assert(c7.get_lhaid(QCDOrder::NLO)  == LhaID(305, 4422, 1, 0));
        assert(c7.get_lhaid(QCDOrder::NNLO) == LhaID(305, 4422, 2, 0));
    }

    {
        SpyProxy sp;
        sp.ret = 7.0;
        const auto v = c7.get_matching_value("LO", ContributionType::SM, std::make_shared<SpyProxy>(sp));
        assert(std::abs(v.real() - 7.0) < 1e-12);
        assert(std::abs(v.imag()) < 1e-12);

        auto spy = std::make_shared<SpyProxy>();
        (void)c7.get_matching_value("NLO", ContributionType::SM, spy);
        assert(spy->last_block == c7.get_storage_block());
        assert(spy->last_id    == c7.id(QCDOrder::NLO, ContributionType::SM));
    }

    std::cout << " INTEGRATION OK\n";
    return 0;
}

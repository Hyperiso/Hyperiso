#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include "SMParamSetter.h"
#include "IMartyParameterProxy.h"
#include "Include.h"

class DummySMProxy : public IMartyParameterProxy<std::string, LhaID> {
public:
    scalar_t operator()(const std::string& blk, const LhaID& id) const override {
        if (blk=="MASS" && id==LhaID(13)) return 0.105;     // m_mu
        if (blk=="MASS" && id==LhaID(3))  return 0.005;     // m_s
        if (blk=="MASS_EW_SCALE" && id==LhaID(5,1)) return 4.7;     // mb(muW)
        if (blk=="MASS_EW_SCALE" && id==LhaID(6))   return 173.0;   // mt(muW)
        if (blk=="SMINPUTS" && id==LhaID(7,1))      return 0.231;   // sin^2(thetaW)

        if (blk=="SMC" && id==LhaID(2)) return scalar_t(1.0, 1.0);
        return 0.0;
    }
};

class DummyBSMProxy : public IMartyParameterProxy<std::string, LhaID> {
public:
    scalar_t operator()(const std::string& blk, const LhaID& id) const override {
        if (blk=="MINPAR" && id==LhaID(3)) return 10.0;        // tanbeta
        if (blk=="XBLK"   && id==LhaID(1)) return scalar_t(3.0,4.0); // complex
        return 0.0;
    }
};

static InterpretedParam P(std::string b, LhaID c, bool bsm, bool cpx){
    return InterpretedParam{std::move(b), c, bsm, cpx};
}

int main(){
    std::cout << "== SMParamSetter UNIT ==\n";

    auto sm  = std::make_shared<DummySMProxy>();
    auto bsm = std::make_shared<DummyBSMProxy>();
    std::set<std::string> specials = {"KIN","WEIN","REGPROP","BETA"};
    SMParamSetter setter("THDM", specials, sm, bsm);

    // KIN 34 → -m_mu^2
    {
        auto m = setter.setParam("kin_mu", P("KIN", LhaID(34), false, false));
        assert(m.count("kin_mu"));
        double exp = -std::pow(0.105,2);
        assert(std::abs(m["kin_mu"] - exp) < 1e-9);
    }
    // WEIN → asin(sqrt(SMINPUTS(7,1)))
    {
        auto m = setter.setParam("wein", P("WEIN", LhaID(7,1), false, false));
        double exp = std::asin(std::sqrt(0.231));
        assert(std::abs(m["wein"] - exp) < 1e-12);
    }
    // MASS 5 and 6 → MASS_EW_SCALE
    {
        auto m5 = setter.setParam("mb", P("MASS", LhaID(5), false, false));
        assert(std::abs(m5["mb"] - 4.7) < 1e-12);
        auto m6 = setter.setParam("mt", P("MASS", LhaID(6), false, false));
        assert(std::abs(m6["mt"] - 173.0) < 1e-12);
    }
    // BSM complex
    {
        auto mc = setter.setParam("zb", P("XBLK", LhaID(1), true, true));
        assert(std::abs(mc["zb_rel"] - 3.0) < 1e-12);
        assert(std::abs(mc["zb_img"] - 4.0) < 1e-12);
    }
    // SM complex
    {
        auto mc = setter.setParam("smc", P("SMC", LhaID(2), false, true));
        assert(std::abs(mc["smc_rel"] - 1.0) < 1e-12);
        assert(std::abs(mc["smc_img"] - 1.0) < 1e-12);
    }

    std::cout << "UNIT OK\n";
    return 0;
}

#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <memory>
#include <set>
#include "SMParamSetter.h"
#include "IMartyParameterProxy.h"
#include "Include.h"

class DummySMProxy : public IMartyParameterProxy<std::string, LhaID> {
public:
    scalar_t operator()(const std::string& blk, const LhaID& id) const override {
        if (blk=="MASS" && id==LhaID(13)) return 0.105;     // m_mu
        if (blk=="MASS" && id==LhaID(1))  return 0.0047;    // m_d
        if (blk=="MASS" && id==LhaID(2))  return 0.00216;   // m_u
        if (blk=="MASS" && id==LhaID(3))  return 0.005;     // m_s
        if (blk=="MASS" && id==LhaID(4))  return 1.273;     // m_c
        if (blk=="MASS_EW_SCALE" && id==LhaID(5,1)) return 4.7;     // mb(muW)
        if (blk=="MASS_EW_SCALE" && id==LhaID(6))   return 173.0;   // mt(muW)
        if (blk=="SMINPUTS" && id==LhaID(7,1))      return 0.231;   // sin^2(thetaW)
        if (blk=="SMINPUTS" && id==LhaID(1))        return 129.39248302366735; // alpha_em^-1

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
    // 2->2 neutral-meson mixing uses a physical threshold point rather
    // than the historical b->s mu mu fallback.  For
    // b + anti-s -> anti-b + s at threshold, p_i.p_j = m_i*m_j.
    {
        const std::string template_path = "/tmp/hyperiso_param_setter_mixing_kinematics.cpp";
        {
            std::ofstream out(template_path);
            out << R"CPP(
                auto wil = model.computeWilsonCoefficients(
                    mty::Order::OneLoop,
                    {Incoming("b"), Incoming(AntiPart("s")),
                     Outgoing(AntiPart("b")), Outgoing("s")}, opts);
            )CPP";
        }

        SMParamSetter mixing_setter("ZPrime", specials, sm, bsm, template_path);
        const double mb = 4.7;
        const double ms = 0.005;

        auto s12 = mixing_setter.setParam("s_12", P("KIN", LhaID(12), false, false));
        auto s13 = mixing_setter.setParam("s_13", P("KIN", LhaID(13), false, false));
        auto s14 = mixing_setter.setParam("s_14", P("KIN", LhaID(14), false, false));
        auto s23 = mixing_setter.setParam("s_23", P("KIN", LhaID(23), false, false));
        auto s24 = mixing_setter.setParam("s_24", P("KIN", LhaID(24), false, false));
        auto s34 = mixing_setter.setParam("s_34", P("KIN", LhaID(34), false, false));

        assert(std::abs(s12["s_12"] - mb*ms) < 1e-12);
        assert(std::abs(s13["s_13"] - mb*mb) < 1e-12);
        assert(std::abs(s14["s_14"] - mb*ms) < 1e-12);
        assert(std::abs(s23["s_23"] - ms*mb) < 1e-12);
        assert(std::abs(s24["s_24"] - ms*ms) < 1e-12);
        assert(std::abs(s34["s_34"] - mb*ms) < 1e-12);

        // Backward compatibility with the historical sm.json KIN aliases.
        auto s12_legacy = mixing_setter.setParam("s_12_old", P("KIN", LhaID(4), false, false));
        auto s13_legacy = mixing_setter.setParam("s_13_old", P("KIN", LhaID(7), false, false));
        assert(std::abs(s12_legacy["s_12_old"] - mb*ms) < 1e-12);
        assert(std::abs(s13_legacy["s_13_old"] - mb*mb) < 1e-12);

        std::remove(template_path.c_str());
    }

    // WEIN → asin(sqrt(SMINPUTS(7,1)))
    {
        auto m = setter.setParam("wein", P("WEIN", LhaID(7,1), false, false));
        double exp = std::asin(std::sqrt(0.231));
        assert(std::abs(m["wein"] - exp) < 1e-12);
    }
    // GAUGE 4 -> e_em = sqrt(4*pi/alpha_em^-1)
    {
        auto m = setter.setParam("e_em", P("GAUGE", LhaID(4), false, false));
        const double expected = std::sqrt(4.0 * std::acos(-1.0) / 129.39248302366735);
        assert(std::abs(m["e_em"] - expected) < 1e-12);
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

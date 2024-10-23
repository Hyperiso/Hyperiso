#include <iostream>
#include "../../ExternalIntegration/MARTY/MARTY_INSTALL/include/marty/models/sm.h"
#include "../../ExternalIntegration/MARTY/MARTY_INSTALL/include/marty.h"

using namespace csl;
using namespace mty;
using namespace std;
using namespace sm_input;

void defineLibPath(Library &lib) {
#ifdef MARTY_LIBRARY_PATH
    lib.addLPath(MARTY_LIBRARY_PATH);
    lib.addLPath(MARTY_LIBRARY_PATH "/..");
    lib.addLPath(MARTY_LIBRARY_PATH "/marty");
    lib.addLPath(MARTY_LIBRARY_PATH "/marty/lha");
#endif
#ifdef MARTY_INCLUDE_PATH
    lib.addIPath(MARTY_INCLUDE_PATH);
#endif
}

vector<Wilson> getO6(const Model& model, const WilsonSet& wilsons, const std::string& up_quark) {
    auto p = wilsons.kinematics.getOrderedMomenta();
    auto b = model.getParticle("b");
    auto s = model.getParticle("s");
    auto q = model.getParticle(up_quark);
    auto gamma = dirac4.gamma;
    auto gamma5 = dirac4.gamma_chir;
    auto i = model.generateIndices(4, "C");
    auto A = model.generateIndex("C", "G");
    auto al = DiracIndices(9);
    auto mu = MinkowskiIndices(3);
    auto T = model.getGenerator("C", "b");

    Expr O6_vv = GetComplexConjugate(s({i[0], al[0]}, p[1])) * gamma({+mu[0], al[0], al[1]}) * gamma({+mu[1], al[1], al[2]}) * gamma({+mu[2], al[2], al[3]}) * T({A, i[0], i[1]}) * b({i[1], al[3]}, p[0]) * 
                    GetComplexConjugate(q({i[2], al[4]}, p[2])) * gamma({mu[0], al[4], al[5]}) * gamma({mu[1], al[5], al[6]}) * gamma({mu[2], al[6], al[7]}) * T({A, i[2], i[3]}) * q({i[3], al[7]}, p[3]);
    Expr O6_av = GetComplexConjugate(s({i[0], al[0]}, p[1])) * gamma({+mu[0], al[0], al[1]}) * gamma({+mu[1], al[1], al[2]}) * gamma({+mu[2], al[2], al[3]}) * gamma5({al[3], al[8]}) * T({A, i[0], i[1]}) * b({i[1], al[3]}, p[0]) * 
                    GetComplexConjugate(q({i[2], al[4]}, p[2])) * gamma({mu[0], al[4], al[5]}) * gamma({mu[1], al[5], al[6]}) * gamma({mu[2], al[6], al[7]}) * T({A, i[2], i[3]}) * q({i[3], al[7]}, p[3]);

    return {{WilsonCoefficient(CSL_HALF), WilsonOperator(O6_vv)}, {WilsonCoefficient(-CSL_HALF), WilsonOperator(O6_av)}};
}

int calculate_C6(Model &model, gauge::Type gauge) {

    model.getParticle("W")->setGaugeChoice(gauge);
    model.getParticle("Z")->setGaugeChoice(gauge);

    undefineNumericalValues(); // Allow for HIso to set all the parameters' values
    mty::option::excludeExternalLegsCorrections = true;

    Expr factorOperator = -4 * GetComplexConjugate(V_ts) * V_tb * G_F / csl::sqrt_s(2);
    FeynOptions opts;
    opts.setFermionOrder({1, 0, 2, 3});
    opts.setWilsonOperatorCoefficient(factorOperator);

    auto wil_u = model.computeWilsonCoefficients(mty::Order::TreeLevel,
        {Incoming("b"), Outgoing("s"),
         Outgoing("u"), Outgoing(AntiPart("u"))},
        opts);

    auto O6_u = getO6(model, wil_u, "u");
    Expr C6_u = getWilsonCoefficient(wil_u, O6_u);
    Replace(C6_u, e_em, sqrt_s(8 * G_F / sqrt_s(2)) * M_W * sin_s(theta_W));

    auto wil_c = model.computeWilsonCoefficients(mty::Order::TreeLevel,
        {Incoming("b"), Outgoing("s"),
         Outgoing("c"), Outgoing(AntiPart("c"))},
        opts);

    auto O6_c = getO6(model, wil_c, "c");
    Expr C6_c = getWilsonCoefficient(wil_c, O6_c);
    Replace(C6_c, e_em, sqrt_s(8 * G_F / sqrt_s(2)) * M_W * sin_s(theta_W));

    [[maybe_unused]] int sysres = system("rm -rf libs/C6_SM");
    mty::Library wilsonLib("C6_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("C6", C6_u + C6_c);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_C6(sm, gauge::Type::Feynman);
}
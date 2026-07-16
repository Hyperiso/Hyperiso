#include <iostream>

// HYPERISO_MARTY_OPERATOR_NORM_ABI: ew-input-normalization-v1
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

vector<Wilson> getO5(const Model& model, const WilsonSet& wilsons, const std::string& up_quark) {
    auto p = wilsons.kinematics.getOrderedMomenta();
    auto b = model.getParticle("b");
    auto s = model.getParticle("s");
    auto q = model.getParticle(up_quark);
    auto gamma = dirac4.gamma;
    auto gamma5 = dirac4.gamma_chir;
    auto i = model.generateIndices(2, "C", s);
    auto al = DiracIndices(9);
    auto mu = MinkowskiIndices(3);

    Expr O5_vv = GetComplexConjugate(s({i[0], al[0]}, p[1])) * gamma({+mu[0], al[0], al[1]}) * gamma({+mu[1], al[1], al[2]}) * gamma({+mu[2], al[2], al[3]}) * b({i[0], al[3]}, p[0]) * 
                    GetComplexConjugate(q({i[1], al[4]}, p[2])) * gamma({mu[0], al[4], al[5]}) * gamma({mu[1], al[5], al[6]}) * gamma({mu[2], al[6], al[7]}) * q({i[1], al[7]}, p[3]);
    Expr O5_av = GetComplexConjugate(s({i[0], al[0]}, p[1])) * gamma({+mu[0], al[0], al[1]}) * gamma({+mu[1], al[1], al[2]}) * gamma({+mu[2], al[2], al[3]}) * gamma5({al[3], al[8]}) * b({i[0], al[8]}, p[0]) * 
                    GetComplexConjugate(q({i[1], al[4]}, p[2])) * gamma({mu[0], al[4], al[5]}) * gamma({mu[1], al[5], al[6]}) * gamma({mu[2], al[6], al[7]}) * q({i[1], al[7]}, p[3]);

    return {{WilsonCoefficient(CSL_HALF), WilsonOperator(O5_vv)}, {WilsonCoefficient(-CSL_HALF), WilsonOperator(O5_av)}};
}

int calculate_C5(Model &model, gauge::Type gauge) {

    model.getParticle("W")->setGaugeChoice(gauge);
    model.getParticle("Z")->setGaugeChoice(gauge);

    undefineNumericalValues(); // Allow for HIso to set all the parameters' values
    mty::option::excludeExternalLegsCorrections = true;

    Expr factorOperator = -GetComplexConjugate(V_ts) * V_tb * pow_s(e_em, 2)
                          / (2 * pow_s(sin_s(theta_W), 2) * pow_s(M_W, 2));
    FeynOptions opts;
    opts.setFermionOrder({1, 0, 2, 3});
    opts.setWilsonOperatorCoefficient(factorOperator);

    auto wil_u = model.computeWilsonCoefficients(mty::Order::TreeLevel,
        {Incoming("b"), Outgoing("s"),
         Outgoing("u"), Outgoing(AntiPart("u"))},
        opts);

    auto O5_u = getO5(model, wil_u, "u");
    Expr C5_u = getWilsonCoefficient(wil_u, O5_u);

    auto wil_c = model.computeWilsonCoefficients(mty::Order::TreeLevel,
        {Incoming("b"), Outgoing("s"),
         Outgoing("c"), Outgoing(AntiPart("c"))},
        opts);

    auto O5_c = getO5(model, wil_c, "c");
    Expr C5_c = getWilsonCoefficient(wil_c, O5_c);

    [[maybe_unused]] int sysres = system("rm -rf libs/C5_SM");
    mty::Library wilsonLib("C5_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("C5", C5_u + C5_c);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_C5(sm, gauge::Type::Feynman);
}
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

int calculate_C4(Model &model, gauge::Type gauge) {

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

    auto O4_u = dimension6Operator(model, wil_u, DiracCoupling::VL, DiracCoupling::V, {"C", ColorCoupling::Generator}, {0, 2, 1, 3});
    Expr C4_u = getWilsonCoefficient(wil_u, O4_u);

    auto wil_c = model.computeWilsonCoefficients(mty::Order::TreeLevel,
        {Incoming("b"), Outgoing("s"),
         Outgoing("c"), Outgoing(AntiPart("c"))},
        opts);

    auto O4_c = dimension6Operator(model, wil_u, DiracCoupling::VL, DiracCoupling::V, {"C", ColorCoupling::Generator}, {0, 2, 1, 3});
    Expr C4_c = getWilsonCoefficient(wil_c, O4_c);

    [[maybe_unused]] int sysres = system("rm -rf libs/C4_SM");
    mty::Library wilsonLib("C4_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("C4", C4_u + C4_c);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_C4(sm, gauge::Type::Feynman);
}
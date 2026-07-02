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

int calculate_C3(Model &model, gauge::Type gauge) {

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

    auto O3_u = dimension6Operator(model, wil_u, DiracCoupling::VL, DiracCoupling::V, {0, 2, 1, 3});
    Expr C3_u = getWilsonCoefficient(wil_u, O3_u);

    auto wil_c = model.computeWilsonCoefficients(mty::Order::TreeLevel,
        {Incoming("b"), Outgoing("s"),
         Outgoing("c"), Outgoing(AntiPart("c"))},
        opts);

    auto O3_c = dimension6Operator(model, wil_u, DiracCoupling::VL, DiracCoupling::V, {0, 2, 1, 3});
    Expr C3_c = getWilsonCoefficient(wil_c, O3_c);

    [[maybe_unused]] int sysres = system("rm -rf libs/C3_SM");
    mty::Library wilsonLib("C3_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("C3", C3_u + C3_c);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_C3(sm, gauge::Type::Feynman);
}
#include <iostream>
// HYPERISO_MARTY_TEMPLATE_ABI: scalar-bqll-finite-scheme-v12

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

int calculate_CQ1e(Model &model, gauge::Type gauge) {

    model.getParticle("W")->setGaugeChoice(gauge);
    model.getParticle("Z")->setGaugeChoice(gauge);

    undefineNumericalValues(); // Allow for HIso to set all the parameters' values
    mty::option::excludeExternalLegsCorrections = false;

    Expr factorOperator = -4 * GetComplexConjugate(V_ts) * V_tb * G_F * pow_s(e_em / (4 * CSL_PI), 2) / csl::sqrt_s(2);
    FeynOptions opts;
    opts.setFermionOrder({1, 0, 2, 3});
    opts.setWilsonOperatorCoefficient(factorOperator);

    auto wil = model.computeWilsonCoefficients(mty::Order::OneLoop,
        {Incoming("b"), Outgoing("s"),
         Outgoing("e"), Outgoing(AntiPart("e"))},
        opts);

    auto Q1 = dimension6Operator(model, wil, DiracCoupling::R, DiracCoupling::S, {1, 0, 2, 3});
    Expr CQ1_e = getWilsonCoefficient(wil, Q1);

    [[maybe_unused]] int sysres = system("rm -rf libs/CQ1_E_SM");
    mty::Library wilsonLib("CQ1_E_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("CQ1_E", CQ1_e);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_CQ1e(sm, gauge::Type::Feynman);
}
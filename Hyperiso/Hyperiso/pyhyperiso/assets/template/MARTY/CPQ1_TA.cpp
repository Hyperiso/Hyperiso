#include <iostream>
// HYPERISO_MARTY_TEMPLATE_ABI: scalar-bqll-oneloop-externallegs-v1

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

int calculate_CPQ1tau(Model &model, gauge::Type gauge) {

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
         Outgoing("tau"), Outgoing(AntiPart("tau"))},
        opts);

    auto QP1 = dimension6Operator(model, wil, DiracCoupling::L, DiracCoupling::S, {1, 0, 2, 3});
    Expr CQP1_tau = getWilsonCoefficient(wil, QP1);

    [[maybe_unused]] int sysres = system("rm -rf libs/CPQ1_TA_SM");
    mty::Library wilsonLib("CPQ1_TA_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("CPQ1_TA", CQP1_tau);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_CPQ1tau(sm, gauge::Type::Feynman);
}
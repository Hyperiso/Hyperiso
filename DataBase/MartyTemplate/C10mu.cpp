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

int calculate_C10mu(Model &model, gauge::Type gauge) {

    model.getParticle("W")->setGaugeChoice(gauge);
    model.getParticle("Z")->setGaugeChoice(gauge);

    undefineNumericalValues(); // Allow for HIso to set all the parameters' values
    mty::option::excludeExternalLegsCorrections = true;

    Expr factorOperator = -4 * GetComplexConjugate(V_ts) * V_tb * G_F * pow_s(e_em / (4 * CSL_PI), 2) / csl::sqrt_s(2);
    FeynOptions opts;
    opts.setFermionOrder({1, 0, 2, 3});
    opts.setWilsonOperatorCoefficient(factorOperator);

    auto wil = model.computeWilsonCoefficients(mty::Order::OneLoop,
        {Incoming("b"), Outgoing("s"),
         Outgoing("mu"), Outgoing(AntiPart("mu"))},
        opts);

    auto O10_mu = dimension6Operator(model, wil, DiracCoupling::VL, DiracCoupling::A, {0, 2, 1, 3});
    Expr C10_mu = getWilsonCoefficient(wil, O10_mu);

    [[maybe_unused]] int sysres = system("rm -rf libs/C10_mu_SM");
    mty::Library wilsonLib("C10_mu_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("C10_mu", C10_mu);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_C10mu(sm, gauge::Type::Feynman);
}
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

int calculate_C2(Model &model, gauge::Type gauge) {

    model.getParticle("W")->setGaugeChoice(gauge);
    model.getParticle("Z")->setGaugeChoice(gauge);

    undefineNumericalValues(); // Allow for HIso to set all the parameters' values
    mty::option::excludeExternalLegsCorrections = true;

    Expr factorOperator = -4 * GetComplexConjugate(V_cs) * V_cb * G_F / csl::sqrt_s(2);
    FeynOptions opts;
    opts.setFermionOrder({1, 3, 2, 0});
    opts.setWilsonOperatorCoefficient(factorOperator);

    auto wil = model.computeWilsonCoefficients(mty::Order::TreeLevel,
        {Incoming("b"), Outgoing("s"),
         Outgoing("c"), Outgoing(AntiPart("c"))},
        opts);

    auto O2 = dimension6Operator(model, wil, DiracCoupling::VL, DiracCoupling::VL, {1, 3, 0, 2});
    Expr C2 = getWilsonCoefficient(wil, O2);
    Replace(C2, e_em, sqrt_s(8 * G_F / sqrt_s(2)) * M_W * sin_s(theta_W));

    [[maybe_unused]] int sysres = system("rm -rf libs/C2_SM");
    mty::Library wilsonLib("C2_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("C2", C2);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_C2(sm, gauge::Type::Feynman);
}
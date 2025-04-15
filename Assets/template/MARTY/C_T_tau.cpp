#include <iostream>
#include "/home/theo/hyperiso/Third_party/MARTY/MARTY_INSTALL/include/marty/models/sm.h"
#include "/home/theo/hyperiso/Third_party/MARTY/MARTY_INSTALL/include/marty.h"

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

int calculate_C_T_tau(Model &model, gauge::Type gauge) {

    model.getParticle("W")->setGaugeChoice(gauge);
    model.getParticle("Z")->setGaugeChoice(gauge);

    undefineNumericalValues(); // Allow for HIso to set all the parameters' values
    mty::option::excludeExternalLegsCorrections = true;

    Expr factorOperator = -4 * V_cb * G_F / csl::sqrt_s(2);
    FeynOptions opts;
    opts.setFermionOrder({1, 0, 2, 3});
    opts.setWilsonOperatorCoefficient(factorOperator);

    auto wil = model.computeWilsonCoefficients(mty::Order::TreeLevel,
        {Incoming("b"), Outgoing("c"),
         Outgoing("tau"), Outgoing(AntiPart("nu_tau"))},
        opts);

    auto O = dimension6Operator(model, wil, DiracCoupling::TL, DiracCoupling::TL, {0, 2, 1, 3});
    Expr C = getWilsonCoefficient(wil, O);

    [[maybe_unused]] int sysres = system("rm -rf libs/C_T_tau_SM");
    mty::Library wilsonLib("C_T_tau_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("C_T_tau_SM", C);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_C_T_tau(sm, gauge::Type::Feynman);
}
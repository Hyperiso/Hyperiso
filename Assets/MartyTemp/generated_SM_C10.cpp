#include <iostream>
#include "/home/theo/hyperiso/Third_party/MARTY/src/MARTY/src/marty/models/sm.h"
#include "/home/theo/hyperiso/Third_party/MARTY/MARTY_INSTALL/include/marty.h"
//42

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

    auto wil = model.computeWilsonCoefficients(
        mty::Order::OneLoop,
        {Incoming("b"), 
         Outgoing("s"),
         Outgoing("mu"),
         Outgoing(AntiPart("mu"))},
        opts
    );

    Expr C10_mu = getWilsonCoefficient(
        wil, 
        dimension6Operator(model, wil, DiracCoupling::VL, DiracCoupling::A)
    ) / 3;

    [[maybe_unused]] int sysres = system("rm -rf libs/C10_SM");
    mty::Library wilsonLib("C10_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("C10", C10_mu);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_C10mu(sm, gauge::Type::Feynman);
}

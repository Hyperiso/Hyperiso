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

int calculate_C7(Model &model, gauge::Type gauge) {

    model.getParticle("W")->setGaugeChoice(gauge);
    model.getParticle("Z")->setGaugeChoice(gauge);

    undefineNumericalValues(); // Allow for HIso to set all the parameters' values
    mty::option::excludeExternalLegsCorrections = true;

    Expr factorOperator = -4 * G_F * GetComplexConjugate(V_ts) * V_tb * e_em * m_b / (16 * CSL_PI * CSL_PI * csl::sqrt_s(2));
    FeynOptions opts;
    opts.setWilsonOperatorCoefficient(factorOperator);

    auto wil = model.computeWilsonCoefficients(mty::Order::OneLoop, 
        {Incoming("b"), Outgoing("s"), Outgoing("A")}, 
        opts);

    auto O7 = chromoMagneticOperator(model, wil, DiracCoupling::R);
    Expr C7 = getWilsonCoefficient(wil, O7);

    [[maybe_unused]] int sysres = system("rm -rf libs/C7_SM");
    mty::Library wilsonLib("C7_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("C7", C7);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_C7(sm, gauge::Type::Feynman);
}
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

int calculate_CP8(Model &model, gauge::Type gauge) {

    model.getParticle("W")->setGaugeChoice(gauge);
    model.getParticle("Z")->setGaugeChoice(gauge);

    undefineNumericalValues(); // Allow for HIso to set all the parameters' values
    mty::option::excludeExternalLegsCorrections = true;

    Expr factorOperator = -4 * G_F * GetComplexConjugate(V_ts) * V_tb * g_s * m_b / (16 * CSL_PI * CSL_PI * csl::sqrt_s(2));
    FeynOptions opts;
    opts.setWilsonOperatorCoefficient(factorOperator);

    auto wil_t = model.computeWilsonCoefficients(mty::Order::OneLoop, 
        {Incoming("b"), Outgoing("s"), Outgoing("G")}, 
        opts);

    auto OP8 = chromoMagneticOperator(model, wil_t, DiracCoupling::L);
    Expr CP8 = getWilsonCoefficient(wil_t, OP8);

    [[maybe_unused]] int sysres = system("rm -rf libs/CP8_SM");
    mty::Library wilsonLib("CP8_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("CP8", CP8);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_CP8(sm, gauge::Type::Feynman);
}
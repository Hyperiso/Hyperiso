#include <iostream>
// HYPERISO_MARTY_GENERIC_TREE_FIRST_SIGNATURE_ABI: v2

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

int calculate_C7(Model &model, gauge::Type gauge) {

    model.getParticle("W")->setGaugeChoice(gauge);
    model.getParticle("Z")->setGaugeChoice(gauge);

    undefineNumericalValues(); // Allow for HIso to set all the parameters' values
    mty::option::excludeExternalLegsCorrections = true;

    Expr factorOperator = -GetComplexConjugate(V_ts) * V_tb * pow_s(e_em, 3) * m_b
                          / (32 * CSL_PI * CSL_PI * pow_s(sin_s(theta_W), 2) * pow_s(M_W, 2));
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
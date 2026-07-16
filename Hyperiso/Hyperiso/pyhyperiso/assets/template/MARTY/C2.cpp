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

int calculate_C2(Model &model, gauge::Type gauge) {

    model.getParticle("W")->setGaugeChoice(gauge);
    model.getParticle("Z")->setGaugeChoice(gauge);

    undefineNumericalValues(); // Allow for HIso to set all the parameters' values
    mty::option::excludeExternalLegsCorrections = true;

    Expr factorOperator = GetComplexConjugate(V_cs) * V_cb * pow_s(e_em, 2)
                          / (2 * pow_s(sin_s(theta_W), 2) * pow_s(M_W, 2));
    FeynOptions opts;
    opts.setFermionOrder({1, 3, 2, 0});
    opts.setWilsonOperatorCoefficient(factorOperator);
    opts.addFilter(mty::filter::disableParticle("G"));

    auto wil_LO = model.computeWilsonCoefficients(
        mty::Order::TreeLevel,
        {Incoming("b"), Outgoing("s"),
         Outgoing("c"), Outgoing(AntiPart("c"))},
        opts
    );

    Show(wil_LO);
    Display(wil_LO);
    Expr C2_LO = getWilsonCoefficient(wil_LO, dimension6Operator(model, wil_LO, DiracCoupling::VL, DiracCoupling::VL, {1, 3, 0, 2}));

    Display(C2_LO);
    // opts.discardLowerOrders = true;
    // auto wil_NLO = model.computeWilsonCoefficients(
    //     mty::Order::OneLoop,
    //     {Incoming("b"), Outgoing("s"),
    //     Outgoing("c"), Outgoing(AntiPart("c"))},
    //     opts
    // );

    // Expr C2_NLO = getWilsonCoefficient(
    //     wil_NLO,
    //     dimension6Operator(model, wil_NLO, DiracCoupling::VL, DiracCoupling::VL)
    // );
    
    // Expr C2 = C2_LO + C2_NLO;

    Expr C2 = C2_LO;

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
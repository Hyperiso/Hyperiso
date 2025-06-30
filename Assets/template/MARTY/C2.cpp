#include <iostream>

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

    Expr factorOperator = 4 * GetComplexConjugate(V_cs) * V_cb * G_F / csl::sqrt_s(2);
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

    Expr C2_LO = getWilsonCoefficient(wil_LO, dimension6Operator(model, wil_LO, DiracCoupling::VL, DiracCoupling::VL, {1, 3, 0, 2}));

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
#include <iostream>
// HYPERISO_MARTY_GENERIC_TREE_FIRST_SIGNATURE_ABI: v2

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

int calculate_C_BD_5(Model &model, gauge::Type gauge) {
    model.getParticle("W")->setGaugeChoice(gauge);
    model.getParticle("Z")->setGaugeChoice(gauge);

    undefineNumericalValues(); // Allow for HIso to set all the parameters' values
    mty::option::excludeExternalLegsCorrections = true;

    Expr factorOperator = 1;
    FeynOptions opts;
    opts.setWilsonOperatorCoefficient(factorOperator);
    opts.setFermionOrder({1, 0, 3, 2});
    opts.addFilter([&](mty::FeynmanDiagram const &diag) { return !diag.contains("G", mty::FeynmanDiagram::DiagramParticleType::Loop);});

    auto wil_t = model.computeWilsonCoefficients(mty::Order::OneLoop, 
        {Incoming("b"), Incoming(AntiPart("d")), Outgoing(AntiPart("b")), Outgoing("d")}, 
        opts);

    // A colour-singlet neutral vector with both L and R FCNC couplings
    // generates the BMU vector-LR operator Q1^LR.  In the SUSY basis used by
    // HyperIso the coefficient transformation is
    //
    //     C_BMU(Q1^LR) = -1/2 * C_SUSY(Q5),
    //
    // hence C5 = -2 C_BMU(Q1^LR).
    //
    // Project the vector-LR structure directly; projecting the crossed scalar
    // Q5 operator does not expose a tree-level vector exchange in MARTY.
    auto O_LR = dimension6Operator(
        model,
        wil_t,
        mty::DiracCoupling::VL,
        mty::DiracCoupling::VR
    );
    // For the [1,0,3,2] four-identical-flavour ordering used here, MARTY's
    // extracted vector-LR coefficient is one quarter of the standard BMU
    // C1^LR normalization.  Therefore:
    //
    //   C5^SUSY = -2 * C1^LR,BMU
    //           = -8 * C_LR,MARTY.
    Expr C = -8 * getWilsonCoefficient(wil_t, O_LR);

    [[maybe_unused]] int sysres = system("rm -rf libs/C_BD_5_SM");
    mty::Library wilsonLib("C_BD_5_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("C_BD_5", C);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_C_BD_5(sm, gauge::Type::Feynman);
}
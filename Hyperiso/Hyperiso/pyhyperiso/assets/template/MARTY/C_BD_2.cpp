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

int calculate_C_BD_2(Model &model, gauge::Type gauge) {
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

    auto O = dimension6Operator(model, wil_t, mty::DiracCoupling::L, mty::DiracCoupling::L);
    Expr C = getWilsonCoefficient(wil_t, O);

    [[maybe_unused]] int sysres = system("rm -rf libs/C_BD_2_SM");
    mty::Library wilsonLib("C_BD_2_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("C_BD_2", C);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_C_BD_2(sm, gauge::Type::Feynman);
}
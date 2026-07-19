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

int calculate_CT_BD_3(Model &model, gauge::Type gauge) {
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

    auto O = dimension6Operator(model, wil_t, mty::DiracCoupling::R, mty::DiracCoupling::R, {"C", mty::ColorCoupling::Crossed});
    Expr C = getWilsonCoefficient(wil_t, O);

    [[maybe_unused]] int sysres = system("rm -rf libs/CT_BD_3_SM");
    mty::Library wilsonLib("CT_BD_3_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("CT_BD_3", C);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_CT_BD_3(sm, gauge::Type::Feynman);
}
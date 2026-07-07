#include <iostream>

using namespace csl;
using namespace mty;
using namespace std;
using namespace sm_input;


namespace {
bool hyperiso_marty_is_light_up_photon_penguin(mty::FeynmanDiagram const& diag) {
    const bool has_light_up_loop =
        diag.contains("u",   mty::FeynmanDiagram::DiagramParticleType::Loop) ||
        diag.contains("u_L", mty::FeynmanDiagram::DiagramParticleType::Loop) ||
        diag.contains("u_R", mty::FeynmanDiagram::DiagramParticleType::Loop) ||
        diag.contains("c",   mty::FeynmanDiagram::DiagramParticleType::Loop) ||
        diag.contains("c_L", mty::FeynmanDiagram::DiagramParticleType::Loop) ||
        diag.contains("c_R", mty::FeynmanDiagram::DiagramParticleType::Loop);

    const bool has_photon_mediator =
        diag.contains("A", mty::FeynmanDiagram::DiagramParticleType::Mediator);

    return has_light_up_loop && has_photon_mediator;
}
} // namespace

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

int calculate_C9mu(Model &model, gauge::Type gauge) {

    model.getParticle("W")->setGaugeChoice(gauge);
    model.getParticle("Z")->setGaugeChoice(gauge);

    undefineNumericalValues(); // Allow for HIso to set all the parameters' values
    mty::option::excludeExternalLegsCorrections = false;

    Expr factorOperator = -4 * GetComplexConjugate(V_ts) * V_tb * G_F * pow_s(e_em / (4 * CSL_PI), 2) / csl::sqrt_s(2);
    FeynOptions opts;
    opts.addFilter([](mty::FeynmanDiagram const& diag) {
        return !hyperiso_marty_is_light_up_photon_penguin(diag);
    });
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
    
    Expr C9_mu = getWilsonCoefficient(
        wil, 
        dimension6Operator(model, wil, DiracCoupling::VL, DiracCoupling::V)
    ) / 3;

    [[maybe_unused]] int sysres = system("rm -rf libs/C9_SM");
    mty::Library wilsonLib("C9_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("C9", C9_mu);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_C9mu(sm, gauge::Type::Feynman);
}
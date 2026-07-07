#include <iostream>
#include <string>

// HYPERISO_MARTY_OPERATOR_NORM_ABI: ew-input-normalization-v1
// HYPERISO_MARTY_TEMPLATE_ABI: semileptonic-c9-no-photon-penguin-v11
using namespace csl;
using namespace mty;
using namespace std;
using namespace sm_input;

namespace {

bool hyperiso_marty_is_photon_name(std::string const& name) {
    return name == "A" || name == "A;\\gamma" || name.find("\\gamma") != std::string::npos;
}

bool hyperiso_marty_has_photon_mediator(mty::FeynmanDiagram const& diag) {
    for (const auto& particle : diag.getParticles(mty::FeynmanDiagram::DiagramParticleType::Mediator)) {
        if (hyperiso_marty_is_photon_name(std::string(particle->getName()))) {
            return true;
        }
    }
    return false;
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
    // The photon-penguin part of b -> s l l is not used directly as a MARTY
    // four-fermion coefficient: at q^2 -> 0 it carries the regulated photon
    // propagator and must be replaced by a finite EFT matching term.  The SM
    // branch is supplied by HyperIso/SuperIso analytically; BSM photon pieces
    // are expected to be added through WilsonMatchingPatch if needed.
    opts.addFilter([](mty::FeynmanDiagram const& diag) {
        return !hyperiso_marty_has_photon_mediator(diag);
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
    );

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

#include <iostream>
#include <string>

// HYPERISO_MARTY_OPERATOR_NORM_ABI: ew-input-normalization-v1
// HYPERISO_MARTY_TEMPLATE_ABI: semileptonic-cp9-vector-penguin-only-v12
using namespace csl;
using namespace mty;
using namespace std;
using namespace sm_input;

namespace {

bool hyperiso_marty_is_photon_name(std::string const& name) {
    return name == "A" || name == "A;\\gamma" || name.find("\\gamma") != std::string::npos;
}

bool hyperiso_marty_is_scalar_particle(const mty::Particle& particle) {
    return particle->getSpinDimension() == 1;
}

bool hyperiso_marty_has_forbidden_c9_penguin_mediator(mty::FeynmanDiagram const& diag) {
    for (const auto& particle : diag.getParticles(mty::FeynmanDiagram::DiagramParticleType::Mediator)) {
        const std::string name = std::string(particle->getName());

        // Photon penguins are not used directly as four-fermion C9/CP9
        // coefficients; finite matching pieces must be supplied analytically
        // or through WilsonMatchingPatch.
        if (hyperiso_marty_is_photon_name(name)) {
            return true;
        }

        // Scalar / Goldstone / neutral-Higgs penguins do not define the vector
        // semileptonic operators C9/CP9.  In some models MARTY can connect a
        // massless G0 mediator and generate reg_prop-dependent pieces that then
        // leak into the vector projection.  Those belong to scalar operators or
        // gauge/EFT subtractions, not to C9/CP9.  Boxes are unaffected because
        // they have no 4-fermion penguin mediator in this category.
        if (hyperiso_marty_is_scalar_particle(particle)) {
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

int calculate_CP9mu(Model &model, gauge::Type gauge) {

    model.getParticle("W")->setGaugeChoice(gauge);
    model.getParticle("Z")->setGaugeChoice(gauge);

    undefineNumericalValues(); // Allow for HIso to set all the parameters' values
    mty::option::excludeExternalLegsCorrections = false;

    Expr factorOperator = -4 * GetComplexConjugate(V_ts) * V_tb * G_F * pow_s(e_em / (4 * CSL_PI), 2) / csl::sqrt_s(2);
    FeynOptions opts;
    // Same photon-penguin policy as C9.  CP9 has no universal SM-like finite
    // replacement for BSM models; add a WilsonMatchingPatch when the model
    // contains a genuine photon-penguin contribution to CP9.
    opts.addFilter([](mty::FeynmanDiagram const& diag) {
        return !hyperiso_marty_has_forbidden_c9_penguin_mediator(diag);
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

    Expr CP9_mu = getWilsonCoefficient(
        wil, 
        dimension6Operator(model, wil, DiracCoupling::VR, DiracCoupling::V)
    );

    [[maybe_unused]] int sysres = system("rm -rf libs/CP9_SM");
    mty::Library wilsonLib("CP9_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("CP9", CP9_mu);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_CP9mu(sm, gauge::Type::Feynman);
}

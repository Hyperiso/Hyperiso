#include <iostream>
#include <string>

// HYPERISO_MARTY_OPERATOR_NORM_ABI: ew-input-normalization-v1
// HYPERISO_MARTY_TEMPLATE_ABI: semileptonic-cp10-split-regprop-photon-component-v16
using namespace csl;
using namespace mty;
using namespace std;
using namespace sm_input;

namespace {

enum class HyperisoMartyC9LinkerSelection {
    NonPhotonVector,
    PhotonOnly
};

HyperisoMartyC9LinkerSelection hyperiso_marty_c9_linker_selection =
    HyperisoMartyC9LinkerSelection::NonPhotonVector;

void hyperiso_marty_set_c9_linker_selection(HyperisoMartyC9LinkerSelection selection) {
    hyperiso_marty_c9_linker_selection = selection;
}

bool hyperiso_marty_is_photon_name(std::string const& name) {
    return name == "A" || name == "A;\\gamma" || name.find("\\gamma") != std::string::npos;
}

bool hyperiso_marty_is_scalar_particle(const mty::Particle& particle) {
    return particle->getSpinDimension() == 1;
}

bool hyperiso_marty_is_photon_linker_particle(const mty::Particle& particle) {
    return hyperiso_marty_is_photon_name(std::string(particle->getName()));
}

bool hyperiso_marty_is_forbidden_cp10_linker_particle(const mty::Particle& particle) {
    // In the split-reg_prop policy, raw photon linkers are evaluated in a
    // separate *_A function.  The main CP10 function keeps only the non-photon
    // vector part, which is evaluated numerically with reg_prop = 1e-6.
    if (hyperiso_marty_is_photon_linker_particle(particle)) {
        return true;
    }

    // Scalar / Goldstone linkers should not be projected onto vector/axial
    // semileptonic coefficients.  They belong to scalar operators or to a
    // separate gauge/EFT treatment, not to CP10.
    if (hyperiso_marty_is_scalar_particle(particle)) {
        return true;
    }

    return false;
}

template <typename Predicate>
bool hyperiso_marty_has_linker_matching(mty::FeynmanDiagram const& diag, Predicate predicate) {
    // In MARTY's 4-fermion construction the penguin linker X in b -> s X is
    // often still External when filters are applied, and may appear as a
    // Mediator only after connectAmplitudes().  Check both categories.
    for (const auto& particle : diag.getParticles(mty::FeynmanDiagram::DiagramParticleType::External)) {
        if (predicate(particle)) {
            return true;
        }
    }
    for (const auto& particle : diag.getParticles(mty::FeynmanDiagram::DiagramParticleType::Mediator)) {
        if (predicate(particle)) {
            return true;
        }
    }
    return false;
}

bool hyperiso_marty_has_photon_linker(mty::FeynmanDiagram const& diag) {
    return hyperiso_marty_has_linker_matching(diag, [](const mty::Particle& particle) {
        return hyperiso_marty_is_photon_linker_particle(particle);
    });
}

bool hyperiso_marty_has_forbidden_cp10_linker(mty::FeynmanDiagram const& diag) {
    return hyperiso_marty_has_linker_matching(diag, [](const mty::Particle& particle) {
        return hyperiso_marty_is_forbidden_cp10_linker_particle(particle);
    });
}

bool hyperiso_marty_accept_c9_linker(mty::FeynmanDiagram const& diag) {
    switch (hyperiso_marty_c9_linker_selection) {
        case HyperisoMartyC9LinkerSelection::NonPhotonVector:
            return !hyperiso_marty_has_forbidden_cp10_linker(diag);
        case HyperisoMartyC9LinkerSelection::PhotonOnly:
            return hyperiso_marty_has_photon_linker(diag);
    }
    return !hyperiso_marty_has_forbidden_cp10_linker(diag);
}

// Backward-compatible aliases used by older generated C9/CP9 templates.
bool hyperiso_marty_has_forbidden_c9_penguin_mediator(mty::FeynmanDiagram const& diag) {
    return hyperiso_marty_has_forbidden_cp10_linker(diag);
}

bool hyperiso_marty_has_photon_mediator(mty::FeynmanDiagram const& diag) {
    return hyperiso_marty_has_photon_linker(diag);
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

int calculate_CP10mu(Model &model, gauge::Type gauge) {

    model.getParticle("W")->setGaugeChoice(gauge);
    model.getParticle("Z")->setGaugeChoice(gauge);

    undefineNumericalValues(); // Allow for HIso to set all the parameters' values
    mty::option::excludeExternalLegsCorrections = false;

    Expr factorOperator = -4 * GetComplexConjugate(V_ts) * V_tb * G_F * pow_s(e_em / (4 * CSL_PI), 2) / csl::sqrt_s(2);
    FeynOptions opts;
    // CP10 uses the same split-reg_prop policy as C9/CP9 when generated in
    // MARTY split mode:
    //   - NonPhotonVector -> CP10 / CP10_SM, evaluated with reg_prop = 1e-6.
    //   - PhotonOnly      -> CP10_A / CP10_SM_A, evaluated with reg_prop = 1.
    opts.addFilter([](mty::FeynmanDiagram const& diag) {
        return hyperiso_marty_accept_c9_linker(diag);
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

    Display(wil);
    Show(wil);
    Expr CP10_mu = getWilsonCoefficient(
        wil, 
        dimension6Operator(model, wil, DiracCoupling::VR, DiracCoupling::A)
    );
    Display(wil);
    [[maybe_unused]] int sysres = system("rm -rf libs/CP10_SM");
    mty::Library wilsonLib("CP10_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("CP10", CP10_mu);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_CP10mu(sm, gauge::Type::Feynman);
}
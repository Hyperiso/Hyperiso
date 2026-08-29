#include <iostream>
#include <string>
#include <vector>

// HYPERISO_MARTY_OPERATOR_NORM_ABI: ew-input-normalization-v1
// HYPERISO_MARTY_TEMPLATE_ABI: semileptonic-cp9-tree-first-split-regprop-recipe-v21
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

bool hyperiso_marty_tree_level_matching = false;

void hyperiso_marty_set_semileptonic_order(mty::Order order) {
    hyperiso_marty_tree_level_matching = (order == mty::Order::TreeLevel);
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

bool hyperiso_marty_is_light_up_name(std::string const& name) {
    // For primed semileptonic charged-Higgs coefficients the SuperIso THDM
    // matching is already written as the GIM-subtracted top contribution.
    // MARTY's raw 4-fermion computation keeps u/c charged-Higgs penguins and
    // boxes separately; after factoring V_tb V_ts^* these pieces carry CKM
    // phases and do not correspond to the short-distance THDM C'_9 matching
    // convention.  Veto only light up-type loop fermions; do not veto external
    // leptons ("mu") or BSM particles such as H+.
    return name == "u" || name == "c"
        || name.find("u_L") != std::string::npos
        || name.find("u_R") != std::string::npos
        || name.find("c_L") != std::string::npos
        || name.find("c_R") != std::string::npos;
}

bool hyperiso_marty_has_light_up_loop(mty::FeynmanDiagram const& diag) {
    for (const auto& particle : diag.getParticles(mty::FeynmanDiagram::DiagramParticleType::Loop)) {
        if (hyperiso_marty_is_light_up_name(std::string(particle->getName()))) {
            return true;
        }
    }
    return false;
}

bool hyperiso_marty_is_forbidden_c9_linker_particle(const mty::Particle& particle) {
    // Raw photon linkers A -> l+l- carry the regulated propagator used by
    // MARTY's 4-fermion penguin patch.  They are not a finite C9/CP9 matching
    // coefficient.  Finite photon pieces must be supplied analytically or via
    // WilsonMatchingPatch.  This does not veto particles inside loops, only the
    // linker that MARTY connects to the lepton current.
    if (hyperiso_marty_is_photon_linker_particle(particle)) {
        return true;
    }

    // Neutral scalar/Goldstone linkers can be present as External particles of
    // the penguin sub-amplitude before MARTY connects the 4-fermion graph. They
    // should not be projected onto vector C9/CP9; scalar effects belong to CQ*
    // operators.  Again, this does not remove scalar BSM particles in loops or
    // boxes, only a scalar linker attached to the lepton current.
    if (hyperiso_marty_is_scalar_particle(particle)) {
        return true;
    }

    return false;
}

template <typename Predicate>
bool hyperiso_marty_has_linker_matching(mty::FeynmanDiagram const& diag, Predicate predicate) {
    // In MARTY's 4-fermion coefficient construction, the penguin linker X in
    // b -> s X is still an External particle when diagram filters are applied.
    // After connectAmplitudes(), it may appear as a Mediator.  Check both;
    // checking Mediator only misses A/G0 and leaves reg_prop-dependent terms.
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

bool hyperiso_marty_has_forbidden_c9_linker(mty::FeynmanDiagram const& diag) {
    return hyperiso_marty_has_linker_matching(diag, [](const mty::Particle& particle) {
        return hyperiso_marty_is_forbidden_c9_linker_particle(particle);
    });
}

bool hyperiso_marty_accept_c9_linker(mty::FeynmanDiagram const& diag) {
    if (hyperiso_marty_tree_level_matching) {
        return true;
    }
    switch (hyperiso_marty_c9_linker_selection) {
        case HyperisoMartyC9LinkerSelection::NonPhotonVector:
            return !hyperiso_marty_has_forbidden_c9_linker(diag);
        case HyperisoMartyC9LinkerSelection::PhotonOnly:
            return hyperiso_marty_has_photon_linker(diag);
    }
    return !hyperiso_marty_has_forbidden_c9_linker(diag);
}

// Backward-compatible aliases: older generated caches/templates called these
// names. Keeping aliases prevents one half-updated template from failing to
// compile.  The semantically precise helper is hyperiso_marty_has_photon_linker.
bool hyperiso_marty_has_forbidden_c9_penguin_mediator(mty::FeynmanDiagram const& diag) {
    return hyperiso_marty_has_forbidden_c9_linker(diag);
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


namespace {

struct HyperisoMartyTreeProjectionTerm {
    std::string id;
    double weight;
    std::vector<int> fermion_order;
    std::vector<int> operator_order;
    mty::DiracCoupling left_current;
    mty::DiracCoupling right_current;
    bool lepton_first;
};

const std::vector<HyperisoMartyTreeProjectionTerm>& hyperiso_marty_tree_projection_recipe() {
    static const std::vector<HyperisoMartyTreeProjectionTerm> terms = {
HYPERISO_MARTY_TREE_PROJECTION_TERMS
    };
    return terms;
}

Expr hyperiso_marty_project_tree_recipe(
    Model& model,
    const WilsonSet& default_wil,
    const FeynOptions& base_opts,
    mty::DiracCoupling default_left,
    mty::DiracCoupling default_right
) {
    const auto& recipe = hyperiso_marty_tree_projection_recipe();
    if (recipe.empty()) {
        // Backward-compatible default: use the physical direct current projector
        // already encoded by the coefficient template.  Per-coefficient F/O
        // overrides therefore retain their historical meaning when no recipe is
        // configured.
        return getWilsonCoefficient(
            default_wil,
            hyperiso_marty_dimension6_operator(
                mty::Order::TreeLevel,
                model,
                default_wil,
                default_left,
                default_right
            )
        );
    }

    Expr result = CSL_0;
    for (const auto& term : recipe) {
        FeynOptions term_opts = base_opts;
        term_opts.setFermionOrder(term.fermion_order);
        auto amplitude = model.computeAmplitude(
            mty::Order::TreeLevel,
            {Incoming("b"),
             Outgoing("s"),
             Outgoing("mu"),
             Outgoing(AntiPart("mu"))},
            term_opts
        );
        if (amplitude.empty()) {
            continue;
        }

        // MARTY's TreeLevel amplitude construction may mutate FeynOptions.
        // Reapply the requested F before matching, exactly as in the validated
        // HyperIso TreeLevel helper.
        term_opts.setFermionOrder(term.fermion_order);
        auto term_wil = model.getWilsonCoefficients(
            amplitude,
            term_opts,
            mty::DecompositionMode::Matching
        );

        const auto left = term.lepton_first ? term.right_current : term.left_current;
        const auto right = term.lepton_first ? term.left_current : term.right_current;

        const Expr projected = getWilsonCoefficient(
            term_wil,
            hyperiso_marty_dimension6_operator(
                mty::Order::TreeLevel,
                model,
                term_wil,
                left,
                right,
                term.operator_order
            )
        );
        result += term.weight * projected;
    }
    return DeepRefreshed(result);
}

} // namespace

int calculate_CP9mu(Model &model, gauge::Type gauge) {

    model.getParticle("W")->setGaugeChoice(gauge);
    model.getParticle("Z")->setGaugeChoice(gauge);

    undefineNumericalValues(); // Allow for HIso to set all the parameters' values
    mty::option::excludeExternalLegsCorrections = false;

    Expr factorOperator = -4 * GetComplexConjugate(V_ts) * V_tb * G_F * pow_s(e_em / (4 * CSL_PI), 2) / csl::sqrt_s(2);
    FeynOptions opts;
    // The photon-penguin part of b -> s l l is not used directly as the final
    // C9 four-fermion coefficient because it carries MARTY's regulated photon
    // propagator.  In BSM-split mode the same template is evaluated twice:
    //   - NonPhotonVector -> exported as CP9 and evaluated numerically
    //                       with reg_prop = 1e-6.
    //   - PhotonOnly      -> exported separately as CP9_A and evaluated
    //                       numerically with reg_prop = 1.
    // The numeric wrapper writes CP9 = CP9_non-photon + CP9_A.
    opts.addFilter([](mty::FeynmanDiagram const& diag) {
        return hyperiso_marty_accept_c9_linker(diag);
    });
    opts.addFilter([](mty::FeynmanDiagram const& diag) {
        return !hyperiso_marty_has_light_up_loop(diag);
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
    
    Expr CP9_mu = CSL_0;
    if (hyperiso_marty_tree_level_matching) {
        CP9_mu = hyperiso_marty_project_tree_recipe(
            model, wil, opts, DiracCoupling::VR, DiracCoupling::V
        );
    } else {
        CP9_mu = getWilsonCoefficient(
            wil,
            dimension6Operator(model, wil, DiracCoupling::VR, DiracCoupling::V)
        );
    }
    
    [[maybe_unused]] int sysres = system("rm -rf libs/C9_SM");
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

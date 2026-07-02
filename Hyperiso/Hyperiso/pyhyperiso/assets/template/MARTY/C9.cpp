#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

// HYPERISO_MARTY_OPERATOR_NORM_ABI: ew-input-normalization-v1
using namespace csl;
using namespace mty;
using namespace std;
using namespace sm_input;

// HYPERISO_MARTY_TEMPLATE_ABI: semileptonic-marty-example-aligned-v7-no-global-model-filters

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

Expr hyperiso_semileptonic_factor() {
    return -pow_s(e_em, 4) * GetComplexConjugate(V_ts) * V_tb
           / (32 * CSL_PI * CSL_PI * pow_s(M_W, 2) * pow_s(sin_s(theta_W), 2));
}

std::vector<Insertion> hyperiso_bsll_insertions() {
    return {Incoming("b"), Outgoing("s"), Outgoing("mu"), Outgoing(AntiPart("mu"))};
}


std::vector<mty::Particle> hyperiso_marty_zcorr_non_sm_particles(mty::Model& model) {
    mty::SM_Model sm_reference;
    mty::Model::current = &model;

    std::unordered_set<std::string> sm_particle_names;
    for (const auto& particle : sm_reference.getParticles()) {
        sm_particle_names.emplace(std::string(particle->getName()));
    }

    std::vector<mty::Particle> non_sm_particles;
    for (const auto& particle : model.getParticles()) {
        const std::string name = std::string(particle->getName());
        if (sm_particle_names.find(name) == sm_particle_names.end()) {
            non_sm_particles.push_back(particle);
        }
    }
    return non_sm_particles;
}

bool hyperiso_marty_zcorr_force_non_sm_particles(mty::FeynOptions& opts, mty::Model& model) {
    auto non_sm_particles = hyperiso_marty_zcorr_non_sm_particles(model);
    if (non_sm_particles.empty()) {
        return false;
    }
    opts.addFilter(mty::filter::forceParticles(non_sm_particles));
    return true;
}

Expr hyperiso_z_penguin_rescale_v() {
    return 4 * csl::pow_s(csl::sin_s(theta_W), 2) * csl::pow_s(csl::cos_s(theta_W), 2);
}

int calculate_C9mu(Model &model, gauge::Type gauge) {

    model.getParticle("W")->setGaugeChoice(gauge);
    model.getParticle("Z")->setGaugeChoice(gauge);
    undefineNumericalValues();
    mty::option::excludeExternalLegsCorrections = false;

    Expr factorOperator = hyperiso_semileptonic_factor();

    FeynOptions opts;
    opts.setFermionOrder({1, 0, 2, 3});
    opts.setWilsonOperatorCoefficient(factorOperator);

    auto insertions = hyperiso_bsll_insertions();

    auto extract_4f = [&](const WilsonSet& wil) {
        return getWilsonCoefficient(
            wil,
            dimension6Operator(
                model,
                wil,
                DiracCoupling::VL,
                DiracCoupling::V
            )
        ) / 3;
    };

    auto wil_tree = model.computeWilsonCoefficients(
        mty::Order::TreeLevel,
        insertions,
        opts
    );
    Expr C9_tree = extract_4f(wil_tree);

    // Follow MARTY's C9_SM reference test for the photon penguin: force the
    // photon line explicitly and do not restrict to Topology::Mass.
    bool old_keep_first_mass_A = mty::option::keepOnlyFirstMassInLoop;
    mty::option::keepOnlyFirstMassInLoop = false;
    FeynOptions opts_A = opts;
    opts_A.addFilter(mty::filter::forceParticle("A"));
    auto wil_A = model.computeWilsonCoefficients(
        mty::Order::OneLoop,
        insertions,
        opts_A
    );
    mty::option::keepOnlyFirstMassInLoop = old_keep_first_mass_A;
    Expr C9_A = extract_4f(wil_A);

    // Follow MARTY's MSSM C9 example: the photon-penguin piece is
    // extracted with Topology::Triangle only (no Mass topology) and is not
    // multiplied by an empirical factor.  The forced-non-SM version is kept as
    // a diagnostic to compare directly with BSM-only MARTY examples.
    Expr C9_A_force_non_sm = CSL_0;
    {
        bool old_keep_first_mass_A_non_sm = mty::option::keepOnlyFirstMassInLoop;
        mty::option::keepOnlyFirstMassInLoop = false;
        FeynOptions opts_A_non_sm = opts;
        opts_A_non_sm.addFilter(mty::filter::forceParticle("A"));
        if (hyperiso_marty_zcorr_force_non_sm_particles(opts_A_non_sm, model)) {
            auto wil_A_non_sm = model.computeWilsonCoefficients(
                mty::Order::OneLoop,
                insertions,
                opts_A_non_sm
            );
            C9_A_force_non_sm = extract_4f(wil_A_non_sm);
        }
        mty::option::keepOnlyFirstMassInLoop = old_keep_first_mass_A_non_sm;
    }

    Expr C9_A_force_non_sm_mssm_style = C9_A_force_non_sm;
    Expr C9_A_corrected = C9_A;

    bool old_keep_first_mass_Z = mty::option::keepOnlyFirstMassInLoop;
    mty::option::keepOnlyFirstMassInLoop = true;
    FeynOptions opts_Z = opts;
    opts_Z.setTopology(Topology::Triangle | Topology::Mass);
    opts_Z.addFilter(mty::filter::disableParticle("A"));
    auto wil_Z = model.computeWilsonCoefficients(
        mty::Order::OneLoop,
        insertions,
        opts_Z
    );
    mty::option::keepOnlyFirstMassInLoop = old_keep_first_mass_Z;
    Expr C9_Z = extract_4f(wil_Z);

    Expr C9_Z_force_non_sm = CSL_0;
    {
        bool old_keep_first_mass_Z_non_sm = mty::option::keepOnlyFirstMassInLoop;
        mty::option::keepOnlyFirstMassInLoop = true;
        FeynOptions opts_Z_non_sm = opts;
        opts_Z_non_sm.setTopology(Topology::Triangle | Topology::Mass);
        opts_Z_non_sm.addFilter(mty::filter::disableParticle("A"));
        if (hyperiso_marty_zcorr_force_non_sm_particles(opts_Z_non_sm, model)) {
            auto wil_Z_non_sm = model.computeWilsonCoefficients(
                mty::Order::OneLoop,
                insertions,
                opts_Z_non_sm
            );
            C9_Z_force_non_sm = extract_4f(wil_Z_non_sm);
        }
        mty::option::keepOnlyFirstMassInLoop = old_keep_first_mass_Z_non_sm;
    }

    // The isolated non-SM Z-penguin from MARTY's four-fermion projection carries
    // explicit Z-lepton normalization factors relative to Hyperiso/SuperIso
    // Wilson-coefficient conventions.  Rescale only the forced non-SM part; the
    // SM-like part is left unchanged.  This keeps generic models safe while
    // fixing BSM loop contributions such as the THDM charged-Higgs Z-penguin.
    Expr C9_Z_force_non_sm_hyperiso_norm =
        hyperiso_z_penguin_rescale_v() * C9_Z_force_non_sm;
    Expr C9_Z_corrected =
        C9_Z + (C9_Z_force_non_sm_hyperiso_norm - C9_Z_force_non_sm);

    bool old_keep_first_mass_box = mty::option::keepOnlyFirstMassInLoop;
    mty::option::keepOnlyFirstMassInLoop = true;
    FeynOptions opts_box = opts;
    opts_box.setTopology(Topology::Box);
    opts_box.addFilter(mty::filter::disableParticle("A"));
    opts_box.addFilter(mty::filter::disableParticle("Z"));
    auto wil_box = model.computeWilsonCoefficients(
        mty::Order::OneLoop,
        insertions,
        opts_box
    );
    mty::option::keepOnlyFirstMassInLoop = old_keep_first_mass_box;
    Expr C9_box = extract_4f(wil_box);

    Expr C9_oneloop_corrected = C9_A_corrected + C9_Z_corrected + C9_box;

    // Generic semileptonic strategy inspired by MARTY's MSSM C9/C10 example:
    // tree-level BSM models (Z', leptoquarks, ...) use the tree-level coefficient
    // directly; otherwise use a one-loop four-fermion decomposition with
    //   A-penguin: Triangle only,
    //   Z-penguin: Triangle | Mass with Hyperiso/SuperIso Z-lepton normalization,
    //   box: Box with A/Z disabled.
    Expr C9_mu = (C9_tree != CSL_0) ? C9_tree : C9_oneloop_corrected;

    [[maybe_unused]] int sysres = system("rm -rf libs/C9_SM");

    mty::Library wilsonLib("C9_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("C9", C9_mu);
    wilsonLib.addFunction("C9_tree", C9_tree);
    wilsonLib.addFunction("C9_A_4f", C9_A);
    wilsonLib.addFunction("C9_A_4f_force_A_marty_sm_style", C9_A);
    wilsonLib.addFunction("C9_A_4f_force_non_sm", C9_A_force_non_sm);
    wilsonLib.addFunction("C9_A_4f_force_non_sm_mssm_style", C9_A_force_non_sm_mssm_style);
    wilsonLib.addFunction("C9_A_4f_corrected", C9_A_corrected);
    wilsonLib.addFunction("C9_Z_4f", C9_Z);
    wilsonLib.addFunction("C9_Z_4f_force_non_sm", C9_Z_force_non_sm);
    wilsonLib.addFunction("C9_Z_4f_force_non_sm_hyperiso_norm", C9_Z_force_non_sm_hyperiso_norm);
    wilsonLib.addFunction("C9_Z_4f_corrected", C9_Z_corrected);
    wilsonLib.addFunction("C9_box_4f", C9_box);
    wilsonLib.addFunction("C9_oneloop_corrected", C9_oneloop_corrected);
    wilsonLib.addFunction("C9_oneloop_mssm_style", C9_oneloop_corrected);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_C9mu(sm, gauge::Type::Feynman);
}

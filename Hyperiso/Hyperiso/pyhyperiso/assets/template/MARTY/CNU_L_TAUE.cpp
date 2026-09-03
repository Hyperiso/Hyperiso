#include <iostream>
#include <string>
#include <vector>
// HYPERISO_MARTY_GENERIC_TREE_FIRST_SIGNATURE_ABI: v2
// HYPERISO_MARTY_OPERATOR_NORM_ABI: bnunu-2301.06990-v1

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

namespace {

struct HyperisoMartyBNuNuTreeProjectionTerm {
    std::string id;
    double weight;
    std::vector<int> fermion_order;
    std::vector<int> operator_order;
    mty::DiracCoupling left_current;
    mty::DiracCoupling right_current;
    bool lepton_first;
};

const std::vector<HyperisoMartyBNuNuTreeProjectionTerm>&
hyperiso_marty_bnunu_tree_projection_recipe() {
    static const std::vector<HyperisoMartyBNuNuTreeProjectionTerm> terms = {
HYPERISO_MARTY_TREE_PROJECTION_TERMS
    };
    return terms;
}

Expr hyperiso_marty_project_bnunu_tree_recipe(
    Model& model,
    const WilsonSet& default_wil,
    const FeynOptions& base_opts,
    const std::vector<Insertion>& insertions,
    mty::DiracCoupling default_left,
    mty::DiracCoupling default_right,
    const std::vector<int>& default_operator_order
) {
    const auto& recipe = hyperiso_marty_bnunu_tree_projection_recipe();
    if (recipe.empty()) {
        return getWilsonCoefficient(
            default_wil,
            hyperiso_marty_dimension6_operator(
                mty::Order::TreeLevel,
                model,
                default_wil,
                default_left,
                default_right,
                default_operator_order
            )
        );
    }

    Expr result = CSL_0;
    for (const auto& term : recipe) {
        FeynOptions term_opts = base_opts;
        term_opts.setFermionOrder(term.fermion_order);
        auto amplitude = model.computeAmplitude(
            mty::Order::TreeLevel, insertions, term_opts
        );
        if (amplitude.empty()) {
            continue;
        }

        // computeAmplitude() can alter the fermion pairing.  Restore the
        // requested F before the matching decomposition, as in the validated
        // direct-Z' C9/C10 recipe path.
        term_opts.setFermionOrder(term.fermion_order);
        auto term_wil = model.getWilsonCoefficients(
            amplitude, term_opts, mty::DecompositionMode::Matching
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

int calculate_CNU_L_TAUE(Model &model, gauge::Type gauge) {
    model.getParticle("W")->setGaugeChoice(gauge);
    model.getParticle("Z")->setGaugeChoice(gauge);

    undefineNumericalValues();
    mty::option::excludeExternalLegsCorrections = true;

    // HyperIso follows arXiv:2301.06990:
    //   L_eff = 4 GF/sqrt(2) lambda_t C O,
    //   O_L/R = e^2/(4pi)^2 (sbar gamma P_L/R b)
    //           (nubar gamma (1-gamma5) nu).
    // MARTY's VL projector is gamma P_L, hence the neutrino current in O is
    // 2*VL.  The factor below therefore carries an extra factor two compared
    // with the C9 vector-current normalisation.  Its overall sign follows the
    // existing HyperIso MARTY Hamiltonian convention.  With the charged-lepton
    // O9/O10 conventions this yields C_L=(C9-C10)/2 when the neutrino and
    // charged-lepton left-handed couplings belong to the same SU(2)_L doublet.
    Expr factorOperator = -8 * GetComplexConjugate(V_ts) * V_tb * G_F
                        * pow_s(e_em / (4 * CSL_PI), 2) / csl::sqrt_s(2);

    FeynOptions opts;
    // Direct neutral-current topology, validated for the Z' s-channel.
    opts.setFermionOrder({0, 1, 2, 3});
    opts.setWilsonOperatorCoefficient(factorOperator);

    // Native call is OneLoop so HyperIso's generic tree-first wrapper can
    // probe TreeLevel first and fall back to OneLoop according to the runtime
    // MartyOrderPolicy.  In BSM mode the wrapper also enforces the non-SM
    // diagram-particle filter before either calculation.
    const std::vector<Insertion> hyperiso_bnunu_insertions = {
        Incoming("b"), Outgoing("s"),
        Outgoing("nu_tau"), Outgoing(AntiPart("nu_e"))
    };
    auto wil = model.computeWilsonCoefficients(
        mty::Order::OneLoop, hyperiso_bnunu_insertions, opts);

    Expr C = CSL_0;
    if (hyperiso_marty_order == mty::Order::TreeLevel) {
        C = hyperiso_marty_project_bnunu_tree_recipe(
            model, wil, opts, hyperiso_bnunu_insertions,
            mty::DiracCoupling::VL, mty::DiracCoupling::VL,
            {1, 2, 0, 3}
        );
    } else {
        auto O = dimension6Operator(
            model, wil, mty::DiracCoupling::VL, mty::DiracCoupling::VL,
            {1, 2, 0, 3});
        C = getWilsonCoefficient(wil, O);
    }

    [[maybe_unused]] int sysres = system("rm -rf libs/CNU_L_TAUE_SM");
    mty::Library wilsonLib("CNU_L_TAUE_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("CNU_L_TAUE", C);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_CNU_L_TAUE(sm, gauge::Type::Feynman);
}

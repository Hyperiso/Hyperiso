#include <iostream>
#include <string>
#include <vector>
// HYPERISO_MARTY_GENERIC_TREE_FIRST_SIGNATURE_ABI: v2
// HYPERISO_MARTY_OPERATOR_NORM_ABI: knunu-diagonal-v1

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

struct HyperisoMartyKNuNuTreeProjectionTerm {
    std::string id;
    double weight;
    std::vector<int> fermion_order;
    std::vector<int> operator_order;
    mty::DiracCoupling left_current;
    mty::DiracCoupling right_current;
    bool lepton_first;
};

const std::vector<HyperisoMartyKNuNuTreeProjectionTerm>&
hyperiso_marty_knunu_tree_projection_recipe() {
    static const std::vector<HyperisoMartyKNuNuTreeProjectionTerm> terms = {
HYPERISO_MARTY_TREE_PROJECTION_TERMS
    };
    return terms;
}

Expr hyperiso_marty_project_knunu_tree_recipe(
    Model& model,
    const WilsonSet& default_wil,
    const FeynOptions& base_opts,
    const std::vector<Insertion>& insertions,
    mty::DiracCoupling default_left,
    mty::DiracCoupling default_right,
    const std::vector<int>& default_operator_order
) {
    const auto& recipe = hyperiso_marty_knunu_tree_projection_recipe();
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

int calculate_CKNU_L_MUMU(Model &model, gauge::Type gauge) {
    model.getParticle("W")->setGaugeChoice(gauge);
    model.getParticle("Z")->setGaugeChoice(gauge);

    undefineNumericalValues();
    mty::option::excludeExternalLegsCorrections = true;

    // Rare-kaon convention used by KPinunuDecay:
    //   H_eff = +4 GF/sqrt(2) lambda_t alpha_e/(4pi) C O,
    //   lambda_t = V_ts^* V_td,
    //   O_L/R = (sbar gamma P_L/R d)(nubar gamma (1-gamma5) nu).
    // MARTY's VL projector contains P_L, so the neutrino current contributes
    // the same factor two already accounted for in the BNuNu templates.
    Expr factorOperator = +8 * GetComplexConjugate(V_ts) * V_td * G_F
                        * pow_s(e_em / (4 * CSL_PI), 2) / csl::sqrt_s(2);

    FeynOptions opts;
    // Direct neutral-current topology for s -> d nu anti-nu matching.
    opts.setFermionOrder({0, 1, 2, 3});
    opts.setWilsonOperatorCoefficient(factorOperator);

    // Native call is OneLoop so HyperIso's generic tree-first wrapper can
    // probe TreeLevel first and fall back to OneLoop according to the runtime
    // MartyOrderPolicy.  In BSM mode the wrapper also enforces the non-SM
    // diagram-particle filter before either calculation.
    const std::vector<Insertion> hyperiso_knunu_insertions = {
        Incoming("d"), Outgoing("s"),
        Outgoing("nu_mu"), Outgoing(AntiPart("nu_mu"))
    };
    auto wil = model.computeWilsonCoefficients(
        mty::Order::OneLoop, hyperiso_knunu_insertions, opts);

    Expr C = CSL_0;
    if (hyperiso_marty_order == mty::Order::TreeLevel) {
        C = hyperiso_marty_project_knunu_tree_recipe(
            model, wil, opts, hyperiso_knunu_insertions,
            mty::DiracCoupling::VL, mty::DiracCoupling::VL,
            {1, 2, 0, 3}
        );
    } else {
        auto O = dimension6Operator(
            model, wil, mty::DiracCoupling::VL, mty::DiracCoupling::VL,
            {1, 2, 0, 3});
        C = getWilsonCoefficient(wil, O);
    }

    [[maybe_unused]] int sysres = system("rm -rf libs/CKNU_L_MUMU_SM");
    mty::Library wilsonLib("CKNU_L_MUMU_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("CKNU_L_MUMU", C);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_CKNU_L_MUMU(sm, gauge::Type::Feynman);
}

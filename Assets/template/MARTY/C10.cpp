#include <iostream>

// HYPERISO_MARTY_OPERATOR_NORM_ABI: ew-input-normalization-v1
using namespace csl;
using namespace mty;
using namespace std;
using namespace sm_input;

// HYPERISO_MARTY_TEMPLATE_ABI: semileptonic-c10-tree-first-full-4f-v9

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

Expr hyperiso_c10_operator_factor() {
    // Keep the historical Hyperiso/SuperIso C10 normalization.  In particular,
    // do not use the EW-input factor introduced for the current-current and
    // dipole templates: the MARTY C9/C10 examples normalize the full bsll
    // four-fermion amplitude with this G_F convention before extracting the
    // axial lepton coefficient.
    return -4 * GetComplexConjugate(V_ts) * V_tb * G_F
           * csl::pow_s(e_em / (4 * CSL_PI), 2) / csl::sqrt_s(2);
}

std::vector<Insertion> hyperiso_bsll_insertions() {
    return {Incoming("b"), Outgoing("s"), Outgoing("mu"), Outgoing(AntiPart("mu"))};
}

int calculate_C10mu(Model &model, gauge::Type gauge) {

    model.getParticle("W")->setGaugeChoice(gauge);
    model.getParticle("Z")->setGaugeChoice(gauge);

    undefineNumericalValues();

    // Important for C10: keep MARTY's external-leg corrections in the amplitude.
    // The earlier decomposed template used excludeExternalLegsCorrections=true
    // and then compensated boxes by hand.  This template intentionally goes back
    // to the full four-fermion extraction: no /3 projection factor, no
    // C10_full - 2*C10_box sign surgery.
    mty::option::excludeExternalLegsCorrections = false;

    Expr factorOperator = hyperiso_c10_operator_factor();

    FeynOptions opts;
    opts.setFermionOrder({1, 0, 2, 3});
    opts.setWilsonOperatorCoefficient(factorOperator);

    auto insertions = hyperiso_bsll_insertions();

    auto extract_C10 = [&](const WilsonSet& wil) {
        return getWilsonCoefficient(
            wil,
            dimension6Operator(
                model,
                wil,
                DiracCoupling::VL,
                DiracCoupling::A
            )
        );
    };

    auto wil_tree = model.computeWilsonCoefficients(
        mty::Order::TreeLevel,
        insertions,
        opts
    );
    Expr C10_tree = extract_C10(wil_tree);

    // Tree-first policy: when the selected operator already has a non-zero
    // tree-level coefficient (for example a Z' or leptoquark mediator), do not
    // evaluate OneLoop at all.  This avoids mixing loop corrections into a
    // matching problem whose leading BSM contribution is tree level.
    Expr C10_full = CSL_0;
    if (C10_tree == CSL_0) {
        auto wil_full = model.computeWilsonCoefficients(
            mty::Order::OneLoop,
            insertions,
            opts
        );
        C10_full = extract_C10(wil_full);
    }

    Expr C10_mu = (C10_tree != CSL_0) ? C10_tree : C10_full;
    std::cout << "[MARTY C10] selected order="
              << ((C10_tree != CSL_0) ? "TreeLevel" : "OneLoop")
              << std::endl;

    [[maybe_unused]] int sysres = system("rm -rf libs/C10_SM");

    mty::Library wilsonLib("C10_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("C10", C10_mu);
    wilsonLib.addFunction("C10_tree", C10_tree);
    wilsonLib.addFunction("C10_full", C10_full);
    wilsonLib.addFunction("C10_oneloop_full", C10_full);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_C10mu(sm, gauge::Type::Feynman);
}

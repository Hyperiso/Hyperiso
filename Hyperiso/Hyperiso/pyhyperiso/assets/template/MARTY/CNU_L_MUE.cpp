#include <iostream>
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

int calculate_CNU_L_MUE(Model &model, gauge::Type gauge) {
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
    // existing HyperIso MARTY Hamiltonian convention (and gives C_L=C9 for a
    // pure SU(2)_L-singlet left-left Z' current).
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
    auto wil = model.computeWilsonCoefficients(
        mty::Order::OneLoop,
        {Incoming("b"), Outgoing("s"),
          Outgoing("nu_mu"), Outgoing(AntiPart("nu_e"))},
        opts);

    auto O = dimension6Operator(
        model, wil, mty::DiracCoupling::VL, mty::DiracCoupling::VL,
        {1, 2, 0, 3});
    Expr C = getWilsonCoefficient(wil, O);

    [[maybe_unused]] int sysres = system("rm -rf libs/CNU_L_MUE_SM");
    mty::Library wilsonLib("CNU_L_MUE_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("CNU_L_MUE", C);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_CNU_L_MUE(sm, gauge::Type::Feynman);
}

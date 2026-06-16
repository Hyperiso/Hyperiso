#include <iostream>

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

int calculate_C10mu(Model &model, gauge::Type gauge) {

    model.getParticle("W")->setGaugeChoice(gauge);
    model.getParticle("Z")->setGaugeChoice(gauge);

    undefineNumericalValues();

    // Convention qui matchait le Z-penguin SuperIso
    mty::option::excludeExternalLegsCorrections = true;
    mty::option::verboseAmplitude = true;

    Expr factorOperator =
        -4 * GetComplexConjugate(V_ts) * V_tb * G_F
        * pow_s(e_em / (4 * CSL_PI), 2)
        / csl::sqrt_s(2);

    FeynOptions opts;
    opts.setFermionOrder({1, 0, 2, 3});
    opts.setWilsonOperatorCoefficient(factorOperator);

    auto insertions = std::vector<Insertion>{
        Incoming("b"),
        Outgoing("s"),
        Outgoing("mu"),
        Outgoing(AntiPart("mu"))
    };

    auto extract_C10 = [&](const WilsonSet& wil) {
        return getWilsonCoefficient(
            wil,
            dimension6Operator(
                model,
                wil,
                DiracCoupling::VL,
                DiracCoupling::A
            )
        ) / 3;
    };

    auto wil_full = model.computeWilsonCoefficients(
        mty::Order::OneLoop,
        insertions,
        opts
    );

    Expr C10_full = extract_C10(wil_full);

    FeynOptions opts_box = opts;
    opts_box.setTopology(Topology::Box);

    auto wil_box = model.computeWilsonCoefficients(
        mty::Order::OneLoop,
        insertions,
        opts_box
    );

    Expr C10_box = extract_C10(wil_box);

    // C10_full = C10_penguin + C10_box_wrong
    // C10_fixed = C10_penguin - C10_box_wrong
    //            = C10_full - 2 * C10_box_wrong
    Expr C10_mu = C10_full - 2 * C10_box;

    [[maybe_unused]] int sysres = system("rm -rf libs/C10_SM");

    mty::Library wilsonLib("C10_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("C10", C10_mu);
    defineLibPath(wilsonLib);
    wilsonLib.print();

    return 0;
}

int main() {
    SM_Model sm;
    return calculate_C10mu(sm, gauge::Type::Feynman);
}
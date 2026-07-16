#include <algorithm>
#include "../../MARTY/MARTY_INSTALL/include/marty/models/sm.h"
#include "../../MARTY/MARTY_INSTALL/include/marty.h"


using namespace csl;
using namespace mty;
using namespace std;
using namespace sm_input;

void defineLibPath(mty::Library &lib)
{
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

int calculate(Model &model, gauge::Type gauge)
{
    using namespace mty::sm_input;
    // undefineNumericalValues();
    model.getParticle("W")->setGaugeChoice(gauge);

    mty::FeynOptions options;
    auto res = model.computeAmplitude(OneLoop, {Incoming("b"), Outgoing("s"), Outgoing("A")}, options);

    // Disable NLO QCD contributions
    res = res.filterOut([&](mty::FeynmanDiagram const &diagram) {return !diagram.contains(model.getParticle("G"), mty::FeynmanDiagram::DiagramParticleType::Loop);});

    Expr V_ts_star      = csl::GetComplexConjugate(V_ts);
    Expr factorOperator = -V_ts_star * V_tb * G_F * e_em / (4 * csl::sqrt_s(2) * CSL_PI * CSL_PI);

    options.setWilsonOperatorCoefficient(factorOperator);
    auto wilsonC7 = model.getWilsonCoefficients(res, options);
    Expr CC7 = getWilsonCoefficient(wilsonC7, chromoMagneticOperator(model, wilsonC7, DiracCoupling::R));
    Expr CC7p = getWilsonCoefficient(wilsonC7, chromoMagneticOperator(model, wilsonC7, DiracCoupling::L));

    [[maybe_unused]] int sysres = system("rm -rf libs/C7_SM");
    mty::Library         wilsonLib("C7_SM", "libs");
    wilsonLib.cleanExistingSources();
    wilsonLib.addFunction("C7", m_b * CC7 + m_s * CC7p);
    defineLibPath(wilsonLib);
    wilsonLib.print();
    return 1;
}

int main()
{

    // mty::sm_input::redefineNumericalValues(); // for compatibility
    SM_Model sm;
    // sm.computeFeynmanRules();
    // Display(sm.getFeynmanRules());

    // calculate(sm, gauge::Type::Unitary);
    // calculate(sm, gauge::Type::Lorenz);
    return calculate(sm, gauge::Type::Feynman);
}
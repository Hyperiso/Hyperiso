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

std::vector<size_t> pos;

int calculate(Model &model, gauge::Type gauge)
{
    using namespace mty::sm_input;
    // undefineNumericalValues();
    model.getParticle("W")->setGaugeChoice(gauge);

    mty::FeynOptions options;
    auto res = model.computeAmplitude(Order::OneLoop,
                                 {Incoming("b"), Outgoing("s"), Outgoing("G")},
                                 options);

    // Disable NLO QCD contributions
    res = res.filterOut([&](mty::FeynmanDiagram const &diagram) {return diagram.contains(model.getParticle("G"), mty::FeynmanDiagram::DiagramParticleType::Loop);});

    Expr V_ts_star      = csl::GetComplexConjugate(V_ts);
    Expr factorOperator = -V_ts_star * V_tb * G_F * g_s / (4 * csl::sqrt_s(2) * CSL_PI * CSL_PI);

    options.setWilsonOperatorCoefficient(factorOperator);
    auto wilsonC8 = model.getWilsonCoefficients(res, options);

    Expr CC8 = getWilsonCoefficient(
        wilsonC8, chromoMagneticOperator(model, wilsonC8, DiracCoupling::R));
    Expr CC8p = getWilsonCoefficient(
        wilsonC8, chromoMagneticOperator(model, wilsonC8, DiracCoupling::L));

    [[maybe_unused]] int sysres = system("rm -rf libs/C8_SM");
    csl::LibraryGenerator::setQuadruplePrecision(false);
    mty::Library wilsonLib("C8_SM", "libs");
    wilsonLib.addFunction("C8", m_b * CC8 + m_s * CC8p);
    defineLibPath(wilsonLib);
    wilsonLib.print();
    return 1;
}

int main()
{
    // mty::sm_input::redefineNumericalValues(); // for compatibility
    SM_Model sm;
    return calculate(sm, gauge::Type::Unitary);
}
// This file is part of MARTY.
//
// MARTY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MARTY is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MARTY. If not, see <https://www.gnu.org/licenses/>.

#include "marty.h"
#include "marty/models/sm.h"

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

auto calculate(Model &model, gauge::Type gauge)
{
    // SM_Model model;
    undefineNumericalValues(); 

    // Keep for testing
    model.getParticle("h")->setEnabledInDiagrams(false);
    model.getParticle("u")->setEnabledInDiagrams(false);
    model.getParticle("c")->setEnabledInDiagrams(false);
    mty::option::keepOnlyFirstMassInLoop = false;

    FeynOptions options;
    options.setFermionOrder({1, 0, 2, 3});
    Expr V_ts_star      = csl::GetComplexConjugate(V_ts);
    Expr factorOperator = -csl::pow_s(e_em, 2) * G_F * V_ts_star * V_tb / (4 * sqrt(2) * CSL_PI * CSL_PI);
    options.setWilsonOperatorCoefficient(factorOperator);

    auto res = model.computeWilsonCoefficients(OneLoop,
                                                    {Incoming("b"),
                                                     Outgoing("s"),
                                                     Outgoing("mu"),
                                                     Outgoing(AntiPart("mu"))},
                                                    options);

    Show(res);

    csl::Expr C9 = getWilsonCoefficient(
        res,
        dimension6Operator(model, res, DiracCoupling::VL, DiracCoupling::V));

    // Library::setQuadruplePrecision(true);
    Library wilsonLib("C9_SM", "libs");
    wilsonLib.addFunction("C9", C9);

    [[maybe_unused]] int sysres = system("rm -rf libs/C9_SM");
    defineLibPath(wilsonLib);
    wilsonLib.print();
}

int main()
{
    SM_Model sm;
    calculate(sm, gauge::Type::Feynman);
    return 0;
}
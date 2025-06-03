#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <cmath>
#include "func_c2_thdm.h"
#include "common.h"

#include "params.h"
#include "libcomplexop.h"
#include "cparams.h"
#include "clib_c2_thdm.h"
#include <complex.h>

namespace c2_thdm {

complex_t C2(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.G_F = param.G_F;
    cparam.M_W = param.M_W;
    cparam.m_b = param.m_b;
    cparam.m_c = param.m_c;
    cparam.V_cb = param.V_cb;
    cparam.e_em = param.e_em;
    cparam.s_13 = param.s_13;
    cparam.theta_W = param.theta_W;
    cparam.reg_prop = param.reg_prop;
    cparam.V_cs = param.V_cs.get().real() + _mty_I*param.V_cs.get().imag();
    auto res = c_C2(&cparam);
    return {res.real, res.imag};
}

} // End of namespace c2_thdm

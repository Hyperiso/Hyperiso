#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <cmath>
#include "marty/core/looptools_interface.h"
#include "func_c2_sm.h"
#include "common.h"

#include "params.h"
#include "libcomplexop.h"
#include "cparams.h"
#include "clib_c2_sm.h"
#include <complex.h>

namespace c2_sm {

complex_t C2(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.G_F = param.G_F;
    cparam.M_W = param.M_W;
    cparam.M_Z = param.M_Z;
    cparam.m_b = param.m_b;
    cparam.m_c = param.m_c;
    cparam.m_h = param.m_h;
    cparam.m_s = param.m_s;
    cparam.V_cb = param.V_cb;
    cparam.s_13 = param.s_13;
    cparam.s_24 = param.s_24;
    cparam.Finite = param.Finite;
    cparam.theta_W = param.theta_W;
    cparam.reg_prop = param.reg_prop;
    cparam.V_cs = param.V_cs.get().real() + _mty_I*param.V_cs.get().imag();
    auto res = c_C2(&cparam);
    return {res.real, res.imag};
}

} // End of namespace c2_sm

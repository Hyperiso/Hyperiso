#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <cmath>
#include "marty/core/looptools_interface.h"
#include "func_c7_thdm.h"
#include "common.h"

#include "params.h"
#include "libcomplexop.h"
#include "cparams.h"
#include "clib_c7_thdm.h"
#include <complex.h>

namespace c7_thdm {

complex_t C7(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.G_F = param.G_F;
    cparam.M_W = param.M_W;
    cparam.m_b = param.m_b;
    cparam.m_c = param.m_c;
    cparam.m_s = param.m_s;
    cparam.m_t = param.m_t;
    cparam.m_u = param.m_u;
    cparam.V_cb = param.V_cb;
    cparam.V_tb = param.V_tb;
    cparam.V_us = param.V_us;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_Hp = param.m_Hp;
    cparam.s_12 = param.s_12;
    cparam.theta_W = param.theta_W;
    cparam.V_cs = param.V_cs.get().real() + _mty_I*param.V_cs.get().imag();
    cparam.V_ts = param.V_ts.get().real() + _mty_I*param.V_ts.get().imag();
    cparam.V_ub = param.V_ub.get().real() + _mty_I*param.V_ub.get().imag();
    auto res = c_C7(&cparam);
    return {res.real, res.imag};
}

} // End of namespace c7_thdm

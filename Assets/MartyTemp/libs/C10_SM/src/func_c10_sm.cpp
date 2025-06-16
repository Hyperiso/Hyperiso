#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <cmath>
#include "marty/core/looptools_interface.h"
#include "func_c10_sm.h"
#include "common.h"

#include "params.h"
#include "libcomplexop.h"
#include "cparams.h"
#include "clib_c10_sm.h"
#include <complex.h>

namespace c10_sm {

complex_t C10(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.G_F = param.G_F;
    cparam.M_W = param.M_W;
    cparam.M_Z = param.M_Z;
    cparam.m_b = param.m_b;
    cparam.m_c = param.m_c;
    cparam.m_s = param.m_s;
    cparam.m_t = param.m_t;
    cparam.m_u = param.m_u;
    cparam.V_cb = param.V_cb;
    cparam.V_tb = param.V_tb;
    cparam.V_us = param.V_us;
    cparam.e_em = param.e_em;
    cparam.m_mu = param.m_mu;
    cparam.s_34 = param.s_34;
    cparam.Finite = param.Finite;
    cparam.theta_W = param.theta_W;
    cparam.reg_prop = param.reg_prop;
    cparam.V_cs = param.V_cs.get().real() + _mty_I*param.V_cs.get().imag();
    cparam.V_ts = param.V_ts.get().real() + _mty_I*param.V_ts.get().imag();
    cparam.V_ub = param.V_ub.get().real() + _mty_I*param.V_ub.get().imag();
    auto res = c_C10(&cparam);
    return {res.real, res.imag};
}

} // End of namespace c10_sm

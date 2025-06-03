#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "c_C2.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_C2(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t G_F = param->G_F;
    const creal_t M_W = param->M_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t V_cb = param->V_cb;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t theta_W = param->theta_W;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t IT_0000 = pow(G_F, -1);
    const ccomplex_t IT_0001 = pow(V_cb, -1);
    const ccomplex_t IT_0002 = cpow(conj(V_cs), -1);
    const ccomplex_t IT_0003 = pow(M_W, 2);
    const ccomplex_t IT_0004 = pow(m_b, 2);
    const ccomplex_t IT_0005 = pow(m_c, 2);
    const ccomplex_t IT_0006 = cpow(2*s_13 + IT_0003 + -IT_0004 + -IT_0005 + 
      -reg_prop, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0008;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0008;
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = IT_0006*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    return create_ccomplex_return((0 + _Complex_I*0.353553390593274)*IT_0000
      *IT_0001*IT_0002*IT_0015);
}


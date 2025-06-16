#ifndef CSL_LIB_CPARAM_H_INCLUDED
#define CSL_LIB_CPARAM_H_INCLUDED

#include "ccommon.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct cparam_s {

    creal_t G_F;
    creal_t M_W;
    creal_t M_Z;
    creal_t m_b;
    creal_t m_c;
    creal_t m_s;
    creal_t m_t;
    creal_t m_u;
    creal_t V_cb;
    creal_t V_tb;
    creal_t V_us;
    creal_t e_em;
    creal_t m_mu;
    creal_t s_12;
    creal_t s_34;
    creal_t Finite;
    creal_t theta_W;
    creal_t reg_prop;
    ccomplex_t V_cs;
    ccomplex_t V_ts;
    ccomplex_t V_ub;

} cparam_t;
#ifdef __cplusplus
} // extern "C" block
#endif


#endif

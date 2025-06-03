#ifndef CSL_LIB_CPARAM_H_INCLUDED
#define CSL_LIB_CPARAM_H_INCLUDED

#include "ccommon.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct cparam_s {

    creal_t G_F;
    creal_t M_W;
    creal_t m_b;
    creal_t m_c;
    creal_t V_cb;
    creal_t s_13;
    creal_t reg_prop;
    ccomplex_t V_cs;

} cparam_t;
#ifdef __cplusplus
} // extern "C" block
#endif


#endif

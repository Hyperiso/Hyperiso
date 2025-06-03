#ifndef CSL_LIB_C2_THDM_CCOMMON_H_INCLUDED
#define CSL_LIB_C2_THDM_CCOMMON_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <complex.h>
#include <stdbool.h>

typedef double creal_t;
typedef double _Complex ccomplex_t;

extern const ccomplex_t _mty_I;

typedef struct ccomplex_return_s {
    creal_t real;
    creal_t imag;
} ccomplex_return_t;

ccomplex_return_t create_ccomplex_return(
    ccomplex_t value);

#ifdef __cplusplus
} // extern "C" block
#endif

#endif

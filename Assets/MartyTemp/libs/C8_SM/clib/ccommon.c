#include "ccommon.h"

const ccomplex_t _mty_I = _Complex_I;

ccomplex_return_t create_ccomplex_return(
    ccomplex_t value)
{
    ccomplex_return_t res = {creal(value), cimag(value)};
    return res;
}

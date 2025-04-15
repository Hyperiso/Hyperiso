#ifndef CSL_LIB_C8_SM_G_H_INCLUDED
#define CSL_LIB_C8_SM_G_H_INCLUDED

#include <array>
#include "common.h"
#include "librarytensor.h"
#include "callable.h"
#include "csl/initSanitizer.h"
#include "params.h"
#include "func_c8_sm.h"

namespace c8_sm {


inline std::array<Callable<complex_t, param_t>, 1> f_G = {
    Callable{"C8", C8},
};

inline std::map<std::string, Callable<complex_t, param_t>> fmap_G {
    {"C8", f_G[0]},
};


}
 // End of namespace c8_sm

#endif

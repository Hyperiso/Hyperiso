#ifndef CSL_LIB_C9_THDM_G_H_INCLUDED
#define CSL_LIB_C9_THDM_G_H_INCLUDED

#include <array>
#include "common.h"
#include "librarytensor.h"
#include "callable.h"
#include "csl/initSanitizer.h"
#include "params.h"
#include "func_c9_thdm.h"

namespace c9_thdm {


inline std::array<Callable<complex_t, param_t>, 1> f_G = {
    Callable{"C9", C9},
};

inline std::map<std::string, Callable<complex_t, param_t>> fmap_G {
    {"C9", f_G[0]},
};


}
 // End of namespace c9_thdm

#endif

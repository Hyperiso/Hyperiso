#ifndef CSL_LIB_C10_THDM_G_H_INCLUDED
#define CSL_LIB_C10_THDM_G_H_INCLUDED

#include <array>
#include "common.h"
#include "librarytensor.h"
#include "callable.h"
#include "csl/initSanitizer.h"
#include "params.h"
#include "func_c10_thdm.h"

namespace c10_thdm {


inline std::array<Callable<complex_t, param_t>, 1> f_G = {
    Callable{"C10", C10},
};

inline std::map<std::string, Callable<complex_t, param_t>> fmap_G {
    {"C10", f_G[0]},
};


}
 // End of namespace c10_thdm

#endif

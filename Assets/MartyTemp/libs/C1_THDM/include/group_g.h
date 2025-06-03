#ifndef CSL_LIB_C1_THDM_G_H_INCLUDED
#define CSL_LIB_C1_THDM_G_H_INCLUDED

#include <array>
#include "common.h"
#include "librarytensor.h"
#include "callable.h"
#include "csl/initSanitizer.h"
#include "params.h"
#include "func_c1_thdm.h"

namespace c1_thdm {


inline std::array<Callable<complex_t, param_t>, 1> f_G = {
    Callable{"C1", C1},
};

inline std::map<std::string, Callable<complex_t, param_t>> fmap_G {
    {"C1", f_G[0]},
};


}
 // End of namespace c1_thdm

#endif

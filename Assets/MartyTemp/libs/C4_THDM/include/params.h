#ifndef CSL_LIB_PARAM_H_INCLUDED
#define CSL_LIB_PARAM_H_INCLUDED

#include <map>
#include <array>
#include "common.h"
#include "libcomplexop.h"
#include "csl/initSanitizer.h"

namespace c4_thdm {

struct param_t {

    ///////////////////////////////////////
    // Elementary parameters to be defined 
    ///////////////////////////////////////



    ///////////////////////////////////////
    // Parameters functions of others  
    // through diagonalization or mass 
    // expressions, see updateSpectrum()  
    // in global.h or set them by hand  
    // 
    // And other default parameters  
    ///////////////////////////////////////


    void reset()
    {
        using real_params = std::array<csl::InitSanitizer<real_t>*, 0>;
        using complex_params = std::array<csl::InitSanitizer<complex_t>*, 0>;

        for (auto &par : real_params{})
        {
            par->reset();
        }

        for (auto &par : complex_params{})
        {
            par->reset();
        }
    }

    void print(std::ostream &out = std::cout) const
    {
        using real_params = std::array<csl::InitSanitizer<real_t> const*, 0>;
        using complex_params = std::array<csl::InitSanitizer<complex_t> const*, 0>;

        out << "param_t struct:\n";
        out << "Real parameters\n";
        for (auto const &par : real_params{})
        {
            out << "  -> " << par->getName() << ": ";
            if (par->hasValue()) {
                out << static_cast<double>(par->get()) << '\n';
            }
            else {
                out << "uninitialized" << '\n';
            }
        }

        out << "Complex parameters\n";
        for (auto const &par : complex_params{})
        {
            out << "  -> " << par->getName() << ": ";
            if (par->hasValue()) {
                out << static_cast<double>(MTY_REAL(par->get()));
                out << " + i*" << static_cast<double>(MTY_IMAG(par->get())) << '\n';
            }
            else {
                out << "uninitialized" << '\n';
            }
        }
        out << "\n";
    }

    std::map<std::string, csl::InitSanitizer<real_t>*> realParams {
    };

    std::map<std::string, csl::InitSanitizer<complex_t>*> complexParams {
    };

};


}

#endif

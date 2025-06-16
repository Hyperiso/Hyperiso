#ifndef CSL_LIB_PARAM_H_INCLUDED
#define CSL_LIB_PARAM_H_INCLUDED

#include <map>
#include <array>
#include "common.h"
#include "libcomplexop.h"
#include "csl/initSanitizer.h"

namespace c10_sm {

struct param_t {

    ///////////////////////////////////////
    // Elementary parameters to be defined 
    ///////////////////////////////////////

    csl::InitSanitizer<real_t> G_F { "G_F" };
    csl::InitSanitizer<real_t> M_W { "M_W" };
    csl::InitSanitizer<real_t> M_Z { "M_Z" };
    csl::InitSanitizer<real_t> m_b { "m_b" };
    csl::InitSanitizer<real_t> m_c { "m_c" };
    csl::InitSanitizer<real_t> m_s { "m_s" };
    csl::InitSanitizer<real_t> m_t { "m_t" };
    csl::InitSanitizer<real_t> m_u { "m_u" };
    csl::InitSanitizer<real_t> V_cb { "V_cb" };
    csl::InitSanitizer<real_t> V_tb { "V_tb" };
    csl::InitSanitizer<real_t> V_us { "V_us" };
    csl::InitSanitizer<real_t> e_em { "e_em" };
    csl::InitSanitizer<real_t> m_mu { "m_mu" };
    csl::InitSanitizer<real_t> s_34 { "s_34" };
    csl::InitSanitizer<real_t> Finite { "Finite" };
    csl::InitSanitizer<real_t> theta_W { "theta_W" };
    csl::InitSanitizer<real_t> reg_prop { "reg_prop" };
    csl::InitSanitizer<complex_t> V_cs { "V_cs" };
    csl::InitSanitizer<complex_t> V_ts { "V_ts" };
    csl::InitSanitizer<complex_t> V_ub { "V_ub" };


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
        using real_params = std::array<csl::InitSanitizer<real_t>*, 17>;
        using complex_params = std::array<csl::InitSanitizer<complex_t>*, 3>;

        for (auto &par : real_params{
                &G_F, &M_W, &M_Z, &m_b, &m_c, 
                &m_s, &m_t, &m_u, &V_cb, &V_tb, &V_us, 
                &e_em, &m_mu, &s_34, &Finite, &theta_W, &reg_prop, })
        {
            par->reset();
        }

        for (auto &par : complex_params{&V_cs, &V_ts, &V_ub, })
        {
            par->reset();
        }
    }

    void print(std::ostream &out = std::cout) const
    {
        using real_params = std::array<csl::InitSanitizer<real_t> const*, 17>;
        using complex_params = std::array<csl::InitSanitizer<complex_t> const*, 3>;

        out << "param_t struct:\n";
        out << "Real parameters\n";
        for (auto const &par : real_params{
                &G_F, &M_W, &M_Z, &m_b, &m_c, 
                &m_s, &m_t, &m_u, &V_cb, &V_tb, &V_us, 
                &e_em, &m_mu, &s_34, &Finite, &theta_W, &reg_prop, })
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
        for (auto const &par : complex_params{&V_cs, &V_ts, &V_ub, })
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
        {"G_F", &G_F},
        {"M_W", &M_W},
        {"M_Z", &M_Z},
        {"m_b", &m_b},
        {"m_c", &m_c},
        {"m_s", &m_s},
        {"m_t", &m_t},
        {"m_u", &m_u},
        {"V_cb", &V_cb},
        {"V_tb", &V_tb},
        {"V_us", &V_us},
        {"e_em", &e_em},
        {"m_mu", &m_mu},
        {"s_34", &s_34},
        {"Finite", &Finite},
        {"theta_W", &theta_W},
        {"reg_prop", &reg_prop},
    };

    std::map<std::string, csl::InitSanitizer<complex_t>*> complexParams {
        {"V_cs", &V_cs},
        {"V_ts", &V_ts},
        {"V_ub", &V_ub},
    };

};


}

#endif

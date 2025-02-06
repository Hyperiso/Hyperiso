#ifndef __QCDHELPER_H__
#define __QCDHELPER_H__

#include <array>
#include <string>
#include "Parameters.h"
#include "Math.h"

struct QCDParamCache {
    double alphas_mZ;
    double m_Z;
    double mt_pole;
    double mb_mb;
    std::array<double, 4> light_masses;

    bool cache_valid();
};

struct SpecialMasses {
    double mc_pole;
    double mb_pole;
    double mb_1S;
    double mt_mt;
};

class QCDHelper {
private:
    static inline QCDParamCache param_cache;
    static inline SpecialMasses special_masses;
    static inline std::array<double, 6> lambdas_running;
    static inline double lambda4_mb_pole;
    static inline double lambda6_mt_pole;

    static inline std::string m_b_type {"pole"};
    static inline std::string m_t_type {"pole"};

    static void update_cached_values();
    static double get_lambda(double mu);
    static std::vector<double> getOrderedMasses();
    static double match_lambda(double target_alpha, double Q, int nf);
    static double alpha_s_explicit(double mu, double lambda, int nf);
    static void set_mass_types(std::string m_b_type, std::string m_t_type);
    static double runMass(double mass, double Q_i, double Q_f, int nf);
    static double R(double alpha, int nf);

    static double calc_mc_pole();
    static double calc_mb_pole();
    static double calc_mb_1S();
    static double calc_mt_mt();

    static void Update();

public:
    static constexpr int Nc = 3;
    static constexpr double C_F = (Nc * Nc - 1.) / (2. * Nc);
    static constexpr double C_A = Nc;

    static constexpr std::array<std::array<double, 3>, 6> beta {{{31. / 3, 134. / 3, 2309.8}, 
                                                              {29. / 3, 115. / 3, 1786.7}, 
                                                              {9      , 32      , 1287.6}, 
                                                              {25. / 3, 77. / 3 , 812.7 }, 
                                                              {23. / 3, 58. / 3 , 361.81}, 
                                                              {7      , 13      , -65   }}};
                                                              
    static constexpr std::array<std::array<double, 3>, 6> gamma {{{2, 8.1388, 34.408},
                                                               {2, 7.8611, 29.678}, 
                                                               {2, 7.5833, 24.840}, 
                                                               {2, 7.3055, 19.894}, 
                                                               {2, 7.0277, 14.839}, 
                                                               {2, 6.75  , 9.6773}}};

    static void Init(double alpha_s_mZ, double m_Z, double mt_pole, double mb_mb, double m_c, double m_s, double m_d, double m_u);

    static double alpha_s(double mu, const std::string& mass_b_type = "pole", const std::string& mass_t_type = "pole");
    static double msbar_mass(int pdg_code, double mu, const std::string& mass_b_type = "pole", const std::string& mass_t_type = "pole");

    static int get_nf(double mu);
    static double mass_c_pole();
    static double mass_b_pole();
    static double mass_b_msbar();
    static double mass_b_1S();
    static double mass_t_pole();
    static double mass_t_msbar();
};

#endif // __QCDHELPER_H__

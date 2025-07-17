#ifndef __QCDHELPER_H__
#define __QCDHELPER_H__

#include <array>
#include <string>
#include "Parameters.h"
#include "Math.h"

struct QCDConstants {
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
};

class QCDHelper {
private:
    static inline MassType m_b_type {MassType::POLE};
    static inline MassType m_t_type {MassType::POLE};

    static double get_lambda(double mu, MassType mass_b_type, MassType mass_t_type);
    static std::vector<double> getOrderedMasses(MassType mass_b_type, MassType mass_t_type);
    static double match_lambda(double target_alpha, double Q, int nf);
    static double alpha_s_explicit(double mu, double lambda, int nf);
    static double runMass(double mass, double Q_i, double Q_f, int nf);
    static double R(double alpha, int nf);

    static double calc_mc_pole(double lambda_4);
    static double calc_mb_pole(double lambda_5);
    static double calc_mb_1S(double lambda_4, double mb_pole);
    static double calc_mt_mt(double lambda6_mt_pole, double lambda_5);

public:
    static inline QCDConstants _static_constants{};
    static inline QCDConstants* constants = &_static_constants;

    static void Init();

    static double alpha_s(double mu, MassType mass_b_type = MassType::POLE, MassType mass_t_type = MassType::POLE);
    static double msbar_mass(int pdg_code, double mu, MassType mass_b_type = MassType::POLE, MassType mass_t_type = MassType::POLE);

    static int get_nf(double mu, MassType mass_b_type = MassType::POLE, MassType mass_t_type = MassType::POLE);
};

#endif // __QCDHELPER_H__

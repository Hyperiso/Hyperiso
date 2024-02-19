// EpsilonCalculator.h
#ifndef EPSILONCALCULATOR_H
#define EPSILONCALCULATOR_H

#include <cmath>
#include <vector>

// class Parameters {
// public:
//     double SM, gp, g2, MSOFT_Q, mass_top_pole, mass_b_pole;
//     double A_b, tan_beta, mu_Q, mass_gluino, mass_b1, mass_b2;
//     double inv_alpha_em, M2_Q, mass_t1, mass_t2, A_t, MqL3_Q, MbR_Q, mass_stl, mass_cha1, mass_cha2;
//     std::vector<double> yut, yub, mass_neut;
//     std::vector<std::vector<double>> sbot_mix, charg_Umix, charg_Vmix, stop_mix, neut_mix;
//     std::vector<std::vector<double>> stop_tan_betamix;

//     QCDParameters run;
//     Parameters(); // Constructeur pour initialiser les paramètres
// };

class EpsilonCalculator {
protected:
    Parameters* sm = Parameters::GetInstance();
    Parameters* susy = Parameters::GetInstance(1);
    static EpsilonCalculator* instance;
    EpsilonCalculator();

public:
    explicit EpsilonCalculator();
    
    double epsilon_0();
    double epsilon_2() const;
    double epsilon_b();
    double epsilon_bp();
    double epsilon_0p();
    double epsilon_1p() const;

    EpsilonCalculator(const EpsilonCalculator&) = delete;
    void operator=(const EpsilonCalculator&) = delete;

    // Méthode statique pour accéder à l'instance
    static EpsilonCalculator* GetInstance() {
        if (!EpsilonCalculator::instance) {
            EpsilonCalculator::instance = new EpsilonCalculator();
        }
        return EpsilonCalculator::instance;
    }
};


// class SMChargedHiggsEpsilonCalculator : public EpsilonCalculator {
// public:
//     using EpsilonCalculator::EpsilonCalculator; // Hérite du constructeur de la classe de base

//     double epsilon_0() override;
//     double epsilon_2() const override;
//     double epsilon_b() override;
//     double epsilon_bp() override;
//     double epsilon_0p() override;
//     double epsilon_1p() const override;
// };

#endif // EPSILONCALCULATOR_H

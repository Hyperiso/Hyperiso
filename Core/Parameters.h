#pragma once

#include "QCDParameters.h"
#include "lha_reader.h"

#include <vector>
#include <array>
#include <map>
#include <string>
#include <complex>

typedef std::complex<double> complex_t; 

enum class FlavorParamType {
    LIFETIME,
    DECAY_CONSTANT,
    DECAY_CONSTANT_RATIO
};

class Parameters {
public:
    // double Q {sm.mass_top_pole};

    QCDParameters QCDRunner;
    static Parameters* GetInstance(int index = 0);

    void setScale(double Q);
    double alpha_s(double Q);
    double running_mass(double quarkmass, double Q_init, double Q_end, std::string option_massb = "running", std::string option_masst = "pole");

    double getFlavorParam(FlavorParamType type, const std::string& id);

    double operator()(std::string block, int pdgCode) {
        if (block == "MASS") {
            return masses[pdgCode];
        }
        if (block =="GAUGE") {
            return gauge[pdgCode];
        }
        if (block=="YUKAWA_CH_U") {
            return lambda_u[pdgCode/10][pdgCode%10];
        }
        if (block=="YUKAWA_CH_D") {
            return lambda_d[pdgCode/10][pdgCode%10];
        }
        if (block=="YUKAWA_CH_L") {
            return lambda_l[pdgCode/10][pdgCode%10];
        }
        if (block=="EXTPAR") {
            return extpar[pdgCode];
        }
        if (block=="ALPHA") {
            return alpha;
        }
        if (block=="HMIX") {
            return hmix[pdgCode];
        }
        if (block == "RECKM") {
            return std::real(ckm[pdgCode/10][pdgCode%10]);
        }
        if (block == "IMCKM") {
            return std::imag(ckm[pdgCode/10][pdgCode%10]);
        }
        if (block == "STOPMIX") {
            return stopmix[pdgCode/10][pdgCode%10];
        }
        if (block == "SBOTMIX") {
            return stopmix[pdgCode/10][pdgCode%10];
        }
        if (block == "UMIX") {
            return umix[pdgCode/10][pdgCode%10];
        }
        if (block == "VMIX") {
            return vmix[pdgCode/10][pdgCode%10];
        }
        if (block == "H0MIX") {
            return A0mix[pdgCode/10][pdgCode%10];
        }
        if (block == "A0MIX") {
            return H0mix[pdgCode/10][pdgCode%10];
        }
        if (block == "MSOFT") {
            return msoft[pdgCode];
        }
        if (block == "YU") {
            return yu[pdgCode/10][pdgCode%10];
        }
        if (block == "YD") {
            return yd[pdgCode/10][pdgCode%10];
        }
        if (block == "AU") {
            return au[pdgCode/10][pdgCode%10];
        }
        if (block== "NMIX") {
            return nmix[pdgCode/10][pdgCode%10];
        }
        return 0;
        
    }

private:
    static Parameters* instance[3];
    Parameters(int modelId); // Constructeur pour initialiser les paramètres
    void initSM();
    void initSUSY();
    void initFlavor();

    std::vector<std::vector<double>> lambda_u, lambda_d, lambda_l;
    
    std::map<int, double> minpar;
    std::map<int, double> extpar;

    std::map<int, double> masses;
    std::map<int, double> gauge;
    std::map<int, double> hmix;
    std::map<int, double> msoft;
    double alpha;
    double susy_Q;

    std::array<std::array<double, 2>, 2> stopmix;
    std::array<std::array<double, 2>, 2> staumix;
    std::array<std::array<double, 2>, 2> sbotmix;
    std::array<std::array<double, 2>, 2> umix;
    std::array<std::array<double, 2>, 2> vmix;
    std::array<std::array<double, 4>, 4> nmix;
    std::array<std::array<double, 4>, 4> A0mix;
    std::array<std::array<double, 4>, 4> H0mix;
    std::array<std::array<double, 3>, 3> yu;
    std::array<std::array<double, 3>, 3> yd;
    std::array<std::array<double, 3>, 3> ye;
    std::array<std::array<double, 3>, 3> au;
    std::array<std::array<double, 3>, 3> ad;
    std::array<std::array<double, 3>, 3> ae;
    std::array<std::array<complex_t, 3>, 3> ckm;

    std::map<std::string, double> lifetimes;
    std::map<std::string, double> fconst;

    Parameters(const Parameters&) = delete;
    Parameters& operator=(const Parameters&) = delete;
    Parameters(Parameters&&) noexcept = default;
    Parameters& operator=(Parameters&&) noexcept = default;
  
};

// struct SM {
//     double SM, gp, g2, MSOFT_Q, mass_top_pole, mass_b_pole, mass_b_Q, mass_t_Q;
// };

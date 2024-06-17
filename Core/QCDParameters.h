#pragma once
#include <cmath>
#include <string>
#include <vector>
#include <tuple>
#include "Math.h"
#include "Logger.h"

class QCDParameters {
public:
    
    
    QCDParameters() { Lambda5 = 0.2; }
    QCDParameters(double alpha_Z, double m_Z, double masst_pole, double massb_b, double mass_u, double mass_d, double mass_s, double mass_c);
    QCDParameters& operator=(const QCDParameters& other) {
        if (this != &other) {
            this->mass_t_pole = other.mass_t_pole;
            this->mass_b_pole = other.mass_b_pole;
            this->mass_b_b = other.mass_b_b;
            this->mass_t_t = other.mass_t_t;
            this->mass_u = other.mass_u;
            this->mass_d = other.mass_d;
            this->mass_s = other.mass_s;
            this->mass_c = other.mass_c;
            this->Lambda5 = other.Lambda5;
            this->m_b_type = other.m_b_type;
            this->m_t_type = other.m_t_type;
        }
        return *this;
    }

    
    double runningAlphasCalculation(double Q, std::string option_massb = "pole", std::string option_masst = "pole");
    double running_mass(double quark_mass, double Qinit, double Qfin,  std::string option_massb = "pole", std::string option_masst = "pole");

    double mb_pole();
    double mc_pole();
    double mb_1S();
    double mt_mt();

    double get_mt_mt() {return this->mass_t_t;}
    double get_mt_pole() {return this->mass_t_pole;}
    double get_mb_mb() {return this->mass_b_b;}
    double get_mb_pole() {return this->mass_b_pole;}

private:
    double Lambda5;
    double mass_u;
    double mass_d;
    double mass_s;
    double mass_c;
    double mass_t_pole;
    double mass_b_pole;
    double mass_b_b;
    double mass_t_t;

    std::string m_b_type;
    std::string m_t_type;

    double alphasRunning(double Q, double Lambda, int nf) const;
    double matchLambda(double target_alpha, double Q, int nf);

    void setMassTypes(std::string m_b_type, std::string m_t_type);
    std::vector<double> getOrderedMasses();
    double runMass(double mass, double Q_i, double Q_f, int nf);
    int getNf(double Q);
    std::tuple<double, double, double> getBetas(int nf) const;
    std::tuple<double, double, double> getGammas(int nf) const;
    double R(double alpha, int nf) const;
};

#pragma once
#include <array>
#include "Parameters.h"
#include "QCDHelper.h"

class Wilson_parameters {

    
    static Wilson_parameters* instance;
    Wilson_parameters();
    
    double mu_W;
    double mu;

public:

    std::shared_ptr<Parameters> sm;
    double alphas_muW;
    double alphas_mu;
    double eta_mu;

    double sw2;
    double xt;
    
    double mass_top_muW;
    double mass_b_muW;
    double mass_b_muW_2;
    double mass_c_muW;

    int nf=5;
	double beta0 = 11.-2./3.*nf;
    int gen{2};
    double ml;
	double xt2;
	double xt3;
	double xt4;
	double xh;

    
    std::array<std::array<double, arraySize>, arraySize> U0 = {};
    std::array<std::array<double, arraySize>, arraySize> U1 = {};
    std::array<std::array<double, arraySize>, arraySize> U2 = {};
    std::array<std::array<double, arraySize>, arraySize> V0 = {};
    std::array<std::array<double, arraySize>, arraySize> V1 = {};
    std::array<double, arraySize> etaMuPowers = {};
    std::array<double, arraySize> etaMuPowers2 = {};

    void SetMu(double mu);
    void SetMuW(double mu_W);
    void set_gen(int new_gen) {this->gen = new_gen; this->ml = (*sm)("MASS", 13+2*(this->gen-2));}
    Wilson_parameters(const Wilson_parameters&) = delete;
    void operator=(const Wilson_parameters&) = delete;
    
    static Wilson_parameters* GetInstance();
};


#pragma once
#include <array>
#include "Parameters.h"
class Wilson_parameters {

    
    static Wilson_parameters* instance;
    Wilson_parameters();
    
    double mu_W;
    double mu;

public:
    static constexpr int arraySize {10};

    Parameters* sm = Parameters::GetInstance();
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
    int gen{3};
    double ml;
	double xt2;
	double xt3;
	double xt4;
	double xh;

    std::array<std::array<std::array<double, 10>, 10>, 10> m00;
    std::array<std::array<std::array<double, 10>, 10>, 10> m10;
    std::array<std::array<std::array<double, 10>, 10>, 10> m11;
    std::array<std::array<std::array<double, 10>, 10>, 10> m20;
    std::array<std::array<std::array<double, 10>, 10>, 10> m21;
    std::array<std::array<std::array<double, 10>, 10>, 10> m22;

    std::array<std::array<std::array<double, 10>, 10>, 10> l00;
    std::array<std::array<std::array<double, 10>, 10>, 10> l01;
    std::array<std::array<std::array<double, 10>, 10>, 10> l10;
    std::array<std::array<std::array<double, 10>, 10>, 10> l11;



    std::array<double, arraySize> ai = {14.0 / 23.0, 16.0 / 23.0, 6.0 / 23.0, -12.0 / 23.0, 0.408619, -0.422989, -0.899395, 0.145649, -1.0, -1.0}; 
    std::array<std::array<double, arraySize>, arraySize> U0 = {};
    std::array<std::array<double, arraySize>, arraySize> U1 = {};
    std::array<std::array<double, arraySize>, arraySize> U2 = {};
    std::array<double, arraySize> etaMuPowers = {};
    

    void SetMu(double mu);
    void SetMuW(double mu_W);
    void set_gen(int new_gen) {this->gen = new_gen; this->ml = (*sm)("MASS", 13+2*(this->gen-2));}
    Wilson_parameters(const Wilson_parameters&) = delete;
    void operator=(const Wilson_parameters&) = delete;
    
    static Wilson_parameters* GetInstance();
};


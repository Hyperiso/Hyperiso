#include "QCDParameters.h"

class Wilson_parameters {

    QCDParameters run;
    


    static constexpr int arraySize {11};
    std::array<double, arraySize> C0w = {};
    std::array<double, arraySize> C1w = {};
    std::array<double, arraySize> C2w = {};
    double mu_W;
    std::array<double, arraySize> C0b = {};
    std::array<double, arraySize> C1b = {};
    std::array<double, arraySize> C2b = {};
    double mu;

    double alphas_muW;
    double alphas_mu;


    double C0w7; 
	double C1w7; 
	double C2w7; 

	double C0w8; 
	double C1w8; 
	double C2w8;


    std::array<std::array<std::array<double, 10>, 10>, 10> m00;
    std::array<std::array<std::array<double, 10>, 10>, 10> m10;
    std::array<std::array<std::array<double, 10>, 10>, 10> m11;
    std::array<std::array<std::array<double, 10>, 10>, 10> m20;
    std::array<std::array<std::array<double, 10>, 10>, 10> m21;
    std::array<std::array<std::array<double, 10>, 10>, 10> m22;


    std::array<double, arraySize> ai = {14.0 / 23.0, 16.0 / 23.0, 6.0 / 23.0, -12.0 / 23.0, 0.408619, -0.422989, -0.899395, 0.145649, -1.0, -1.0}; 
    std::array<std::array<double, arraySize>, arraySize> U0 = {};
    std::array<std::array<double, arraySize>, arraySize> U1 = {};
    std::array<std::array<double, arraySize>, arraySize> U2 = {};
    std::array<double, arraySize> etaMuPowers = {};

    public:
    Wilson_parameters();
};
#include "Wilson_susy.h"


void SUSY_LO_Strategy::init(Parameters* sm, double scale, WilsonSet& C_match) {

    // Parameters* sm = Parameters::GetInstance();
    double epsilonbp,epsilon0p,epsilon0,epsilon2,epsilon1p,epsilonb;
	if(param->THDM_model==0)
	{	
		epsilonbp=epsilon_bp(param);
		epsilon0p=epsilon_0p(param);
		epsilon0=epsilon_0(param);
		epsilon2=epsilon_2(param);
		epsilon1p=epsilon_1p(param);
		epsilonb=epsilon0+epsilon2;	
	}
}
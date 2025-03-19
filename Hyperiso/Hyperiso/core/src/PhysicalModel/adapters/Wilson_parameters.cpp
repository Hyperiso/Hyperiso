#include "Wilson_parameters.h"
#include <cmath>

Wilson_parameters::Wilson_parameters() {
	LOG_DEBUG("WilsonParameters creation");
    sm = Parameters::GetInstance(ParameterType::SM);
	alphas_muW=QCDHelper::alpha_s(81);
}


Wilson_parameters* Wilson_parameters::GetInstance() {
        if (!Wilson_parameters::instance) {
            Wilson_parameters::instance = new Wilson_parameters();
        }
        return Wilson_parameters::instance;
    }

void Wilson_parameters::SetMuW(double mu_W) {
	
	this->mu_W = mu_W;
	LOG_DEBUG("mu_W : " + std::to_string(mu_W));
	alphas_muW=QCDHelper::alpha_s(mu_W);
	LOG_DEBUG("ALPHA AFTER CALCULATION :", alphas_muW);
	mass_top_muW=QCDHelper::msbar_mass(6, mu_W, "running"); //mass top at top ?
	LOG_DEBUG("mass_top_muW : " + std::to_string(mass_top_muW));
	mass_b_muW=QCDHelper::msbar_mass(5, mu_W, "running"); //mass bottom 6 (at pole)
	LOG_DEBUG("mass_b_muW : " + std::to_string(mass_b_muW));
	mass_b_muW_2=QCDHelper::msbar_mass(5, mu_W);
	mass_c_muW=QCDHelper::msbar_mass(4, mu_W, "pole");
	sw2=pow(sin(atan((*sm)("GAUGE",1)/(*sm)("GAUGE",2))),2.);
	xt = std::pow(mass_top_muW / (*sm)("MASS",24), 2);
	ml = (*sm)("MASS", 13+2*(this->gen-2));
	nf=5;
	beta0 = 11.-2./3.*nf;
	xt2=xt*xt;
	xt3=xt*xt2;
	xt4=xt*xt3;
	xh=pow((*sm)("MASS",25)/(*sm)("MASS",24),2.);

	LOG_DEBUG("Alpha_s at " +std::to_string(mu_W) +" : " + std::to_string(alphas_muW));
	LOG_DEBUG("Mass of top quark at scale: " + std::to_string(mass_top_muW));
    LOG_DEBUG("Mass of bottom quark at scale: " + std::to_string(mass_b_muW));
	LOG_DEBUG("Square of weak mixing angle: " + std::to_string(sw2));
    LOG_DEBUG("Square of top quark mass over W boson mass (xt): " + std::to_string(xt));

}

void Wilson_parameters::SetMu(double mu) {

	std::shared_ptr<Parameters> sm = Parameters::GetInstance();

	this->mu = mu;
	alphas_mu=QCDHelper::alpha_s(mu);	
	eta_mu=alphas_muW/alphas_mu;

	for (int i = 0; i < arraySize; ++i) {
        (etaMuPowers)[i] = std::pow(eta_mu, (ai)[i]);
    }
	for (int i = 0; i < arraySize; ++i) {
        (etaMuPowers2)[i] = std::pow(eta_mu, (ai2)[i]);
    }
	

	LOG_DEBUG("U0,U1, U2 for a scale of " + std::to_string(mu));
	for (int ke = 0; ke < arraySize; ++ke) {
        for (int le = 0; le < arraySize; ++le) {
            (U0)[ke][le] =0;
            (U1)[ke][le] = 0;
            (U2)[ke][le] = 0;
			V0[ke][le] = 0;
			V1[ke][le] = 0;
            for (int ie = 0; ie < arraySize; ++ie) {
                (U0)[ke][le] += (m00)[ke][le][ie] * (etaMuPowers)[ie];
                (U1)[ke][le] += (m10)[ke][le][ie] * (etaMuPowers)[ie] + (m11)[ke][le][ie] * (etaMuPowers)[ie] / eta_mu;
                (U2)[ke][le] += (m20)[ke][le][ie] * (etaMuPowers)[ie] + (m21)[ke][le][ie] * (etaMuPowers)[ie] / eta_mu + (m22)[ke][le][ie] * (etaMuPowers[ie]) / (eta_mu * eta_mu);

				V0[ke][le]=V0[ke][le] + l00[ke][le][ie]*pow(eta_mu,ai[ie]);
		        V1[ke][le]=V1[ke][le] + l10[ke][le][ie]*pow(eta_mu,ai[ie])+l11[ke][le][ie]*pow(eta_mu,ai[ie]-1.);
            }
			
            LOG_DEBUG("U0[" + std::to_string(ke) + "][" + std::to_string(le) + "]: " + std::to_string((U0)[ke][le]));
            LOG_DEBUG("U1[" + std::to_string(ke) + "][" + std::to_string(le) + "]: " + std::to_string((U1)[ke][le]));
            LOG_DEBUG("U2[" + std::to_string(ke) + "][" + std::to_string(le) + "]: " + std::to_string((U2)[ke][le]));
        }
    }

}


Wilson_parameters* Wilson_parameters::instance = nullptr;
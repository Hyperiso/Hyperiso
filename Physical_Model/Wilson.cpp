#include "Wilson.h"
#include <iostream>

// void WilsonInitializer::init() {
    
// }

std::array<complex_t, 11> extractCoefficients(const WilsonSet& C_match, int order) {
    std::array<complex_t, 11> coefficients = {};

    if (!C_match.empty() && C_match[order].size() >= 10) {
        for (size_t i = 1; i <= 10; ++i) {
            coefficients[i] = C_match[order][i - 1];
        }
    }

    return coefficients;
}


WilsonManager* WilsonManager::instance = nullptr;

void SM_LO_Strategy::init(Parameters* sm, double scale, WilsonSet& C_match) {

    double mass_top_muW=(*sm).run.runningAlphasCalculation((*sm)("MASS",6)); //mass top at top ?
	double mass_b_muW=(*sm).run.runningAlphasCalculation((*sm)("MASS",5)); //mass bottom 6 (at pole)
    // or the opposite

    double L=log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24)); // scale -> mu_W
 	double sw2=pow(sin(atan((*sm)("COUPLING",1)/(*sm)("COUPLING",2))),2.); //1 = param-> gp and 2 = param->g2


	double xt= pow(mass_top_muW/(*sm)("MASS",24),2.); // W boson mass (24)
	double yt= pow(mass_top_muW/(*sm)("MASS",25),2.); // param->mass_H (25)

	/* LO */
	
	double C2SM_0 = 1.;
	double C7SM_0 = -0.5*A0t(xt)-23./36.;
	double C8SM_0 = -0.5*F0t(xt)-1./3.;
	double C9SM_0 = (1.-4.*sw2)/sw2*C0t(xt)-B0t(xt)/sw2-D0t(xt) +1./4./sw2+38./27.-4./9.*L;
	double C10SM_0 = (B0t(xt)-C0t(xt))/sw2-1./4./sw2;

	if (C_match.size() < 1) C_match.resize(1);  
    auto& C_LO = C_match[0];
	C_LO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));

	C_LO[static_cast<size_t>(WilsonCoefficient::C2)] = complex_t(C2SM_0, 0);
    C_LO[static_cast<size_t>(WilsonCoefficient::C7)] = complex_t(C7SM_0, 0);
    C_LO[static_cast<size_t>(WilsonCoefficient::C8)] = complex_t(C8SM_0, 0);
    C_LO[static_cast<size_t>(WilsonCoefficient::C9)] = complex_t(C9SM_0, 0);
    C_LO[static_cast<size_t>(WilsonCoefficient::C10)] = complex_t(C10SM_0, 0);



}


void SM_LO_Strategy::set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) {

	Parameters* sm = Parameters::GetInstance();
	auto C_matchs = extractCoefficients(C_match, 0);


	constexpr double pi = 3.141592654;
	double alphas_muW=(*sm).run.runningAlphasCalculation(Q_match);
	double alphas_mu=(*sm).run.runningAlphasCalculation(Q);	
	double eta_mu=alphas_muW/alphas_mu;

	complex_t C7_eff= C_matchs[7]-1./3.*C_matchs[3]-4./9.*C_matchs[4]-20./3.*C_matchs[5]-80./9.*C_matchs[6]; 
	complex_t C8_eff= C_matchs[8]+C_matchs[3]-1./6.*C_matchs[4]+20.*C_matchs[5]-10./3.*C_matchs[6]; 

	Wilson_parameters *W_param = Wilson_parameters::GetInstance();
	for (int i = 0; i < W_param->arraySize; ++i) {
        (W_param->etaMuPowers)[i] = std::pow(eta_mu, (W_param->ai)[i]);
    }

	for (int ke = 0; ke < W_param->arraySize; ++ke) {
        for (int le = 0; le < W_param->arraySize; ++le) {
            for (int ie = 0; ie < W_param->arraySize; ++ie) {
                (W_param->U0)[ke][le] += (W_param->m00)[ke][le][ie] * (W_param->etaMuPowers)[ie];
                (W_param->U1)[ke][le] += (W_param->m10)[ke][le][ie] * (W_param->etaMuPowers)[ie] + (W_param->m11)[ke][le][ie] * (W_param->etaMuPowers)[ie] / eta_mu;
                (W_param->U2)[ke][le] += (W_param->m20)[ke][le][ie] * (W_param->etaMuPowers)[ie] + (W_param->m21)[ke][le][ie] * (W_param->etaMuPowers)[ie] / eta_mu + (W_param->m22)[ke][le][ie] * (W_param->etaMuPowers[ie]) / (eta_mu * eta_mu);
            }
        }
    }

    auto calculateC0b = [&](int ie, int je) {
        return (W_param->U0)[ie-1][je-1] * (je <= 6 ? C_matchs[je] : (je == 7 ? C_matchs[7] : C_matchs[8]));
    };


	if (C.size() < 1) C.resize(1);  
    auto& C_LO = C[0];
	C_LO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));

	for (int ie = 1; ie <= 8; ie++) {
        for (int je = 1; je <= 8; je++) {
            C_LO[ie-1] += calculateC0b(ie, je);
        }
    }

	double fourPiOverAlphasMu = 4.0 * pi / alphas_mu;

    auto updateC0b = [&](int je) {
        return (W_param->U0)[8][je-1] * C_matchs[je];
    };


    for (int je = 1; je <= 8; je++) {
        C_LO[9-1] += fourPiOverAlphasMu * updateC0b(je);

    }

    C_LO[10-1] = C_LO[10-1];

}

void SM_LO_Strategy::set_base2(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) {

	Parameters* sm = Parameters::GetInstance();
	auto C_matchs = extractCoefficients(C_match, 0);

	double alphas_muW=(*sm).run.runningAlphasCalculation(Q_match); //mt pole and mb pole
	double alphas_mu=(*sm).run.runningAlphasCalculation(Q); //mt pole and mb pole
	double eta_mu=alphas_muW/alphas_mu;
	
	complex_t C0w7= C_matchs[7]-1./3.*C_matchs[5]-C_matchs[6]; 

	complex_t C0w8= C_matchs[8]+C_matchs[5];

	std::vector<double> etaMuPowers;

    for (auto exponent : {6./23., -12./23., 0.4086, -0.4230, -0.8994, 0.1456, 16./23., 14./23., 11./23., 29./23.}) {
        etaMuPowers.push_back(std::pow(eta_mu, exponent));
    }

	if (C.size() < 1) C.resize(1);
    auto& C_LO = C[0];
	C_LO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));

	C_LO[0] = (0.5 * etaMuPowers[0] - 0.5 * etaMuPowers[1]) * C_matchs[2];  // C1
    C_LO[1] = (0.5 * etaMuPowers[0] + 0.5 * etaMuPowers[1]) * C_matchs[2];  // C2
    C_LO[2] = (-1./14. * etaMuPowers[0] + 1./6. * etaMuPowers[1] + 0.0509 * etaMuPowers[2] - 0.1403 * etaMuPowers[3] - 0.01126 * etaMuPowers[4] + 0.0054 * etaMuPowers[5]) * C_matchs[2];  // C3
    C_LO[3] = (-1./14. * etaMuPowers[0] - 1./6. * etaMuPowers[1] + 0.0984 * etaMuPowers[2] + 0.1214 * etaMuPowers[3] + 0.0156 * etaMuPowers[4] + 0.0026 * etaMuPowers[5]) * C_matchs[2];  // C4
    C_LO[4] = (-0.0397 * etaMuPowers[2] + 0.0117 * etaMuPowers[3] - 0.0025 * etaMuPowers[4] + 0.0304 * etaMuPowers[5]) * C_matchs[2];  // C5
    C_LO[5] = (0.0335 * etaMuPowers[2] + 0.0239 * etaMuPowers[3] - 0.0462 * etaMuPowers[4] - 0.0112 * etaMuPowers[5]) * C_matchs[2];  // C6
    C_LO[6] = std::pow(eta_mu, 16./23.) * C_matchs[7] + 8./3. * (std::pow(eta_mu, 14./23.) - std::pow(eta_mu, 16./23.)) * C_matchs[8] + C_matchs[2] * (2.2996 * etaMuPowers[7] - 1.0880 * etaMuPowers[6] - 3./7. * etaMuPowers[0] - 1./14. * etaMuPowers[1] - 0.6494 * etaMuPowers[2] - 0.0380 * etaMuPowers[3] - 0.0185 * etaMuPowers[4] - 0.0057 * etaMuPowers[5]);  // C7
    C_LO[7] = std::pow(eta_mu, 14./23.) * C_matchs[8] + C_matchs[2] * (0.8623 * etaMuPowers[7] - 0.9135 * etaMuPowers[2] + 0.0873 * etaMuPowers[3] - 0.0571 * etaMuPowers[4] + 0.0209 * etaMuPowers[5]);  // C8
    C_LO[8] = C_matchs[9] + 4. * pi / alphas_muW * (-4. / 33. * (1. - etaMuPowers[8]) + 8. / 87. * (1. - etaMuPowers[9])) * C_matchs[2];  // C9
    C_LO[9] = C_matchs[10];  // C10



}


void SM_NLO_Strategy::init(Parameters* sm, double scale, WilsonSet& C_match) {

	double mass_top_muW=(*sm).run.runningAlphasCalculation((*sm)("MASS",6)); //mass top at top ?
	double mass_b_muW=(*sm).run.runningAlphasCalculation((*sm)("MASS",5)); //mass bottom 6 (at pole)
    // or the opposite

	double xt= pow(mass_top_muW/(*sm)("MASS",24),2.); // W boson mass (24)
	double yt= pow(mass_top_muW/(*sm)("MASS",25),2.); // param->mass_H (25)

	double L=log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24)); // scale -> mu_W
 	double sw2=pow(sin(atan((*sm)("COUPLING",1)/(*sm)("COUPLING",2))),2.); //1 = param-> gp and 2 = param->g2


    double C1SM_1 = 15.+6.*L;
	double C4SM_1 = E0t(xt)-7./9.+2./3.*L;
	double C7SM_1 = -0.5*A1t(xt,log(scale*scale/mass_top_muW/mass_top_muW))+713./243.+4./81.*L-4./9.*C4SM_1;
	double C8SM_1 = -0.5*F1t(xt,log(scale*scale/mass_top_muW/mass_top_muW))+91./324.-4./27.*L-C4SM_1/6.;
	double C9SM_1 = (1.-4.*sw2)/sw2*C1t(xt,log(scale*scale/mass_top_muW/mass_top_muW))-B1t(xt,log(scale*scale/mass_top_muW/mass_top_muW))/sw2-D1t(xt,log(scale*scale/mass_top_muW/mass_top_muW)) +1./sw2+524./729.-128./243.*pi*pi-16./3.*L-128./81.*L*L;
	double C10SM_1 = (B1t(xt,log(scale*scale/mass_top_muW/mass_top_muW))-C1t(xt,log(scale*scale/mass_top_muW/mass_top_muW)))/sw2-1./sw2;

	if (C_match.size() < 1) C_match.resize(1);  
    auto& C_NLO = C_match[1];
	C_NLO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));

	C_NLO[static_cast<size_t>(WilsonCoefficient::C1)] = complex_t(C1SM_1, 0);
    C_NLO[static_cast<size_t>(WilsonCoefficient::C4)] = complex_t(C4SM_1, 0);
    C_NLO[static_cast<size_t>(WilsonCoefficient::C7)] = complex_t(C7SM_1, 0);
    C_NLO[static_cast<size_t>(WilsonCoefficient::C8)] = complex_t(C8SM_1, 0);
    C_NLO[static_cast<size_t>(WilsonCoefficient::C9)] = complex_t(C9SM_1, 0);
	C_NLO[static_cast<size_t>(WilsonCoefficient::C10)] = complex_t(C10SM_1, 0);
}


void SM_NLO_Strategy::set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) {

	Parameters* sm = Parameters::GetInstance();
	auto C_matchs = extractCoefficients(C_match, 1);
	auto C0_matchs = extractCoefficients(C_match, 0);

	constexpr double pi = 3.141592654;
	double alphas_muW=(*sm).run.runningAlphasCalculation(Q_match);
	double alphas_mu=(*sm).run.runningAlphasCalculation(Q);	
	double eta_mu=alphas_muW/alphas_mu;

	complex_t C7_eff= C_matchs[7]-1./3.*C_matchs[3]-4./9.*C_matchs[4]-20./3.*C_matchs[5]-80./9.*C_matchs[6]; 
	complex_t C8_eff= C_matchs[8]+C_matchs[3]-1./6.*C_matchs[4]+20.*C_matchs[5]-10./3.*C_matchs[6]; 


	complex_t C7_eff_0= C0_matchs[7]-1./3.*C0_matchs[3]-4./9.*C0_matchs[4]-20./3.*C0_matchs[5]-80./9.*C0_matchs[6]; 
	complex_t C8_eff_0= C0_matchs[8]+C0_matchs[3]-1./6.*C0_matchs[4]+20.*C0_matchs[5]-10./3.*C0_matchs[6]; 


	Wilson_parameters *W_param = Wilson_parameters::GetInstance();
	for (int i = 0; i < W_param->arraySize; ++i) {
        (W_param->etaMuPowers)[i] = std::pow(eta_mu, (W_param->ai)[i]);
    }

	for (int ke = 0; ke < W_param->arraySize; ++ke) {
        for (int le = 0; le < W_param->arraySize; ++le) {
            for (int ie = 0; ie < W_param->arraySize; ++ie) {
                (W_param->U0)[ke][le] += (W_param->m00)[ke][le][ie] * (W_param->etaMuPowers)[ie];
                (W_param->U1)[ke][le] += (W_param->m10)[ke][le][ie] * (W_param->etaMuPowers)[ie] + (W_param->m11)[ke][le][ie] * (W_param->etaMuPowers)[ie] / eta_mu;
                (W_param->U2)[ke][le] += (W_param->m20)[ke][le][ie] * (W_param->etaMuPowers)[ie] + (W_param->m21)[ke][le][ie] * (W_param->etaMuPowers)[ie] / eta_mu + (W_param->m22)[ke][le][ie] * (W_param->etaMuPowers[ie]) / (eta_mu * eta_mu);
            }
        }
    }


    auto calculateC1b = [&](int ie, int je) {
        complex_t u0_term = (W_param->U0)[ie-1][je-1] * (je <= 6 ? C_matchs[je] : (je == 7 ? C7_eff : C8_eff));
        complex_t u1_term = (W_param->U1)[ie-1][je-1] * (je <= 6 ? C0_matchs[je] : (je == 7 ? C7_eff_0 : C8_eff_0));
        return eta_mu * (u0_term + u1_term);
    };


	if (C.size() < 1) C.resize(1);  
    auto& C_NLO = C[1];
	C_NLO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));

	for (int ie = 1; ie <= 8; ie++) {
        for (int je = 1; je <= 8; je++) {
            C_NLO[ie-1] += calculateC1b(ie, je);
        }
    }

	double fourPiOverAlphasMu = 4.0 * pi / alphas_mu;


    auto updateC1b = [&](int je) {
        return eta_mu * ((W_param->U0)[8][je-1] * C_matchs[je] + (W_param->U1)[8][je-1] * C0_matchs[je]);
    };


    for (int je = 1; je <= 8; je++) {

        C_NLO[9] += fourPiOverAlphasMu * updateC1b(je);

    }

    C_NLO[9] += fourPiOverAlphasMu * eta_mu * (W_param->U0)[8][8] * C_matchs[8];

    C_NLO[10] = eta_mu * C_matchs[10];

}

void SM_NLO_Strategy::set_base2(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) {

	Parameters* sm = Parameters::GetInstance();
	auto C_matchs_0 = extractCoefficients(C_match, 0);
	auto C_matchs = extractCoefficients(C_match, 1);

	double alphas_muW=(*sm).run.runningAlphasCalculation(Q_match); //mt pole and mb pole
	double alphas_mu=(*sm).run.runningAlphasCalculation(Q); //mt pole and mb pole
	double eta_mu=alphas_muW/alphas_mu;
	
	complex_t C0w7= C_matchs_0[7]-1./3.*C_matchs_0[5]-C_matchs_0[6]; 
	complex_t C1w7= C_matchs[7]-1./3.*C_matchs[5]-C_matchs[6]; 

	complex_t C0w8= C_matchs_0[8]+C_matchs_0[5];
	complex_t C1w8= C_matchs[8]+C_matchs[5]; 
 
	

	std::vector<double> etaMuPowers;

	for (auto exponent : {6./23., -12./23., 0.4086, -0.4230, -0.8994, 0.1456, 39./23., 37./23., 11./23., 29./23.}) {
        etaMuPowers.push_back(std::pow(eta_mu, exponent));
    }

	if (C.size() < 1) C.resize(1);
    auto& C_NLO = C[1];

	C_NLO[0] = (C_matchs_0[2] * 0.8136 + 1.0197 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[0] + (C_matchs_0[2] * 0.7142 + 2.9524 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[1];  // C1
    C_NLO[1] = (C_matchs_0[2] * 0.8136 + 1.0197 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[0] - (C_matchs_0[2] * 0.7142 + 2.9524 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[1];  // C2
    C_NLO[2] = (-0.0766 * C_matchs_0[2] - 0.1457 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[0] + (-0.1455 * C_matchs_0[2] - 0.9841 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[1]
               + (0.1494 * eta_mu * C_matchs[4] - 0.8848 * C_matchs_0[2] + 0.2303 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[2]
               + (-0.3726 * eta_mu * C_matchs[4] + 0.4137 * C_matchs_0[2] + 1.4672 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[3]
               + (0.0738 * eta_mu * C_matchs[4] - 0.0114 * C_matchs_0[2] + 0.0971 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[4]
               + (-0.0173 * eta_mu * C_matchs[4] + 0.1722 * C_matchs_0[2] - 0.0213 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[5];  // C3

	C_NLO[3] = (-0.2353 * C_matchs_0[2] - 0.1457 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[0]
			+ (-0.0397 * C_matchs_0[2] + 0.9841 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[1]
			+ (0.2885 * eta_mu * C_matchs[4] + 0.4920 * C_matchs_0[2] + 0.4447 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[2]
			+ (0.3224 * eta_mu * C_matchs[4] - 0.2758 * C_matchs_0[2] - 1.2696 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[3]
			+ (-0.1025 * eta_mu * C_matchs[4] + 0.0019 * C_matchs_0[2] - 0.1349 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[4]
			+ (-0.0084 * eta_mu * C_matchs[4] - 0.1449 * C_matchs_0[2] - 0.0104 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[5];

	// C5
	C_NLO[4] = 0.0397 * C_matchs_0[2] * etaMuPowers[0] + 0.0926 * C_matchs_0[2] * etaMuPowers[1]
			+ (-0.1163 * eta_mu * C_matchs[4] + 0.7342 * C_matchs_0[2] - 0.1792 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[2]
			+ (0.0310 * eta_mu * C_matchs[4] - 0.1262 * C_matchs_0[2] - 0.1221 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[3]
			+ (0.0162 * eta_mu * C_matchs[4] - 0.1209 * C_matchs_0[2] + 0.0213 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[4]
			+ (-0.0975 * eta_mu * C_matchs[4] - 0.1085 * C_matchs_0[2] - 0.1197 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[5];

	C_NLO[5] = -0.1191 * C_matchs_0[2] * etaMuPowers[0] - 0.2778 * C_matchs_0[2] * etaMuPowers[1]
           + (0.0982 * eta_mu * C_matchs[4] - 0.5544 * C_matchs_0[2] + 0.1513 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[2]
           + (0.0634 * eta_mu * C_matchs[4] + 0.1915 * C_matchs_0[2] - 0.2497 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[3]
           + (0.3026 * eta_mu * C_matchs[4] - 0.2744 * C_matchs_0[2] + 0.3983 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[4]
           + (0.0358 * eta_mu * C_matchs[4] + 0.3568 * C_matchs_0[2] + 0.0440 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[5];


	C_NLO[6] = std::pow(eta_mu, 39./23.) * C1w7 
           + 8./3. * (std::pow(eta_mu, 37./23.) - std::pow(eta_mu, 39./23.)) * C1w8 
           + (297664./14283. * std::pow(eta_mu, 16./23.) - 7164416./357075. * std::pow(eta_mu, 14./23.) + 256868./14283. * std::pow(eta_mu, 37./23.) - 6698884./357075. * std::pow(eta_mu, 39./23.)) * C_matchs_0[8]
           + 37208./4761. * (std::pow(eta_mu, 39./23.) - std::pow(eta_mu, 16./23.)) * C_matchs_0[7]
           + (4661194./816831. * eta_mu * C_matchs[4] - 17.3023 * C_matchs_0[2] + 14.8088 * eta_mu * C_matchs[1] / 15.) * std::pow(eta_mu, 14./23.)
           + (-8516./2217. * eta_mu * C_matchs[4] + 8.5027 * C_matchs_0[2] - 10.8090 * eta_mu * C_matchs[1] / 15.) * std::pow(eta_mu, 16./23.)
           + (4.5508 * C_matchs_0[2] - 0.8740 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[0]
           + (0.7519 * C_matchs_0[2] + 0.4218 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[1]
           + (-1.9043 * eta_mu * C_matchs[4] + 2.0040 * C_matchs_0[2] - 2.9347 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[2]
           + (-0.1008 * eta_mu * C_matchs[4] + 0.7476 * C_matchs_0[2] + 0.3971 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[3]
           + (0.1216 * eta_mu * C_matchs[4] - 0.5385 * C_matchs_0[2] + 0.1600 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[4]
           + (0.0183 * eta_mu * C_matchs[4] + 0.0914 * C_matchs_0[2] + 0.0225 * eta_mu * C_matchs[1] / 15.) * etaMuPowers[5];

	// C9
	C_NLO[8] = eta_mu * (C_matchs[9] + 4. * pi / alphas_muW * (-4. / 33. * (1. - etaMuPowers[8]) + 8. / 87. * (1. - etaMuPowers[9])) * C_matchs[2]);

	// C10
	C_NLO[9] = eta_mu * C_matchs[10];
}


void SM_NNLO_Strategy::init(Parameters* sm, double scale, WilsonSet& C_match) {


	double mass_top_muW=(*sm).run.runningAlphasCalculation((*sm)("MASS",6)); //mass top at top ?
	double mass_b_muW=(*sm).run.runningAlphasCalculation((*sm)("MASS",5)); //mass bottom 6 (at pole)
    // or the opposite

	double xt= pow(mass_top_muW/(*sm)("MASS",24),2.); // W boson mass (24)
	double yt= pow(mass_top_muW/(*sm)("MASS",25),2.); // param->mass_H (25)

	double L=log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24)); // scale -> mu_W
 	double sw2=pow(sin(atan((*sm)("COUPLING",1)/(*sm)("COUPLING",2))),2.); //1 = param-> gp and 2 = param->g2

    double C1SM_2 = -T(xt)+7987./72.+17.*pi*pi/3.+475./6.*L+17.*L*L;
	double C2SM_2 = 127./18.+4./3.*pi*pi+46./3.*L+4.*L*L;
	double C3SM_2 = G1t(xt,log(scale*scale/mass_top_muW/mass_top_muW))-680./243.-20./81.*pi*pi-68./81.*L-20./27.*L*L;
	double C4SM_2 = E1t(xt,log(scale*scale/mass_top_muW/mass_top_muW))+950./243.+10./81.*pi*pi+124./27.*L+10./27.*L*L;
	double C5SM_2 = -G1t(xt,log(scale*scale/mass_top_muW/mass_top_muW))/10.+2./15.*E0t(xt)+68./243.+2./81.*pi*pi+14./81.*L+2./27.*L*L;
	double C6SM_2 = -3./16.*G1t(xt,log(scale*scale/mass_top_muW/mass_top_muW))+E0t(xt)/4.+85./162.+5./108.*pi*pi+35./108.*L+5./36.*L*L;

	double xtW=pow((*sm).run.runningAlphasCalculation((*sm)("MASS",6))/(*sm)("MASS",24),2.); // mass top at mass top
	double xtt=pow((*sm)("MASS",6)/(*sm)("MASS",24),2.); // 24 -> W
	
	double C7SM_2 = (C7t2mt(xtt)+log(scale*scale/mass_top_muW/mass_top_muW)*((-592.*pow(xt,5.)-22.*pow(xt,4.)+12814.*pow(xt,3.)-6376.*xt*xt+512.*xt)/27./pow(xt-1.,5.)*Li2(1.-1./xt)
	+(-26838.*pow(xt,5.)+25938.*pow(xt,4.)+627367.*pow(xt,3.)-331956.*xt*xt+16989.*xt-460.)/729./pow(xt-1.,6.)*log(xt)
	+(34400.*pow(xt,5.)+276644.*pow(xt,4.)-2668324.*pow(xt,3.)+1694437.*xt*xt-323354.*xt+53077.)/2187./pow(xt-1.,5.)
	+log(scale*scale/mass_top_muW/mass_top_muW)*((-63.*pow(xt,5.)+532.*pow(xt,4.)+2089.*pow(xt,3.)-1118.*xt*xt)/9./pow(xt-1.,6.)*log(xt)
	+(1186.*pow(xt,5.)-2705.*pow(xt,4.)-24791.*pow(xt,3.)-16099.*xt*xt+19229.*xt-2740.)/162./pow(xt-1.,5.))) )
	-(C7c2MW(xtW)+13763./2187.*log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24))+814./729.*pow(log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24)),2.));

	double C8SM_2 = (C8t2mt(xtt)+log(scale*scale/mass_top_muW/mass_top_muW)*((-148.*pow(xt,5.)+1052.*pow(xt,4.)-4811.*pow(xt,3.)-3520.*xt*xt-61.*xt)/18./pow(xt-1.,5.)*Li2(1.-1./xt)
	+(-15984.*pow(xt,5.)+152379.*pow(xt,4.)-1358060.*pow(xt,3.)-1201653.*xt*xt-74190.*xt+9188.)/1944./pow(xt-1.,6.)*log(xt)
	+(109669.*pow(xt,5.)-1112675.*pow(xt,4.)+6239377.*pow(xt,3.)+8967623.*xt*xt+768722.*xt-42796.)/11664./pow(xt-1.,5.)
	+log(scale*scale/mass_top_muW/mass_top_muW)*((-139.*pow(xt,4.)-2938.*pow(xt,3.)-2683.*xt*xt)/12./pow(xt-1.,6.)*log(xt)
	+(1295.*pow(xt,5.)-7009.*pow(xt,4.)+29495.*pow(xt,3.)+64513.*xt*xt+17458.*xt-2072.)/216./pow(xt-1.,5.))) )
	-(C8c2MW(xtW)+16607./5832.*log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24))+397./486.*pow(log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24)),2.));
	
	if (C_match.size() < 1) C_match.resize(1);  
    auto& C_NNLO = C_match[2];
	C_NNLO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));

	C_NNLO[static_cast<size_t>(WilsonCoefficient::C1)] = complex_t(C1SM_2, 0);
    C_NNLO[static_cast<size_t>(WilsonCoefficient::C2)] = complex_t(C2SM_2, 0);
    C_NNLO[static_cast<size_t>(WilsonCoefficient::C3)] = complex_t(C3SM_2, 0);
    C_NNLO[static_cast<size_t>(WilsonCoefficient::C4)] = complex_t(C4SM_2, 0);
    C_NNLO[static_cast<size_t>(WilsonCoefficient::C5)] = complex_t(C5SM_2, 0);
	C_NNLO[static_cast<size_t>(WilsonCoefficient::C6)] = complex_t(C6SM_2, 0);
	C_NNLO[static_cast<size_t>(WilsonCoefficient::C7)] = complex_t(C7SM_2, 0);
	C_NNLO[static_cast<size_t>(WilsonCoefficient::C8)] = complex_t(C8SM_2, 0);


}


void SM_NNLO_Strategy::set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) {

	Parameters* sm = Parameters::GetInstance();

	auto C_matchs = extractCoefficients(C_match, 2);
	auto C1_matchs = extractCoefficients(C_match, 1);
	auto C0_matchs = extractCoefficients(C_match, 0);

	constexpr double pi = 3.141592654;
	double alphas_muW=(*sm).run.runningAlphasCalculation(Q_match);
	double alphas_mu=(*sm).run.runningAlphasCalculation(Q);	
	double eta_mu=alphas_muW/alphas_mu;

	complex_t C7_eff= C_matchs[7]-1./3.*C_matchs[3]-4./9.*C_matchs[4]-20./3.*C_matchs[5]-80./9.*C_matchs[6]; 
	complex_t C8_eff= C_matchs[8]+C_matchs[3]-1./6.*C_matchs[4]+20.*C_matchs[5]-10./3.*C_matchs[6]; 

	complex_t C7_eff_0= C0_matchs[7]-1./3.*C0_matchs[3]-4./9.*C0_matchs[4]-20./3.*C0_matchs[5]-80./9.*C0_matchs[6]; 
	complex_t C8_eff_0= C0_matchs[8]+C0_matchs[3]-1./6.*C0_matchs[4]+20.*C0_matchs[5]-10./3.*C0_matchs[6]; 

	complex_t C7_eff_1= C1_matchs[7]-1./3.*C1_matchs[3]-4./9.*C1_matchs[4]-20./3.*C1_matchs[5]-80./9.*C1_matchs[6]; 
	complex_t C8_eff_1= C1_matchs[8]+C1_matchs[3]-1./6.*C1_matchs[4]+20.*C1_matchs[5]-10./3.*C1_matchs[6]; 


	Wilson_parameters *W_param = Wilson_parameters::GetInstance();
	for (int i = 0; i < W_param->arraySize; ++i) {
        (W_param->etaMuPowers)[i] = std::pow(eta_mu, (W_param->ai)[i]);
    }

	for (int ke = 0; ke < W_param->arraySize; ++ke) {
        for (int le = 0; le < W_param->arraySize; ++le) {
            for (int ie = 0; ie < W_param->arraySize; ++ie) {
                (W_param->U0)[ke][le] += (W_param->m00)[ke][le][ie] * (W_param->etaMuPowers)[ie];
                (W_param->U1)[ke][le] += (W_param->m10)[ke][le][ie] * (W_param->etaMuPowers)[ie] + (W_param->m11)[ke][le][ie] * (W_param->etaMuPowers)[ie] / eta_mu;
                (W_param->U2)[ke][le] += (W_param->m20)[ke][le][ie] * (W_param->etaMuPowers)[ie] + (W_param->m21)[ke][le][ie] * (W_param->etaMuPowers)[ie] / eta_mu + (W_param->m22)[ke][le][ie] * (W_param->etaMuPowers[ie]) / (eta_mu * eta_mu);
            }
        }
    }


    auto calculateC2b = [&](int ie, int je) {
        complex_t u0_term = (W_param->U0)[ie-1][je-1] * (je <= 6 ? C_matchs[je] : (je == 7 ? C7_eff : C8_eff));
        complex_t u1_term = (W_param->U1)[ie-1][je-1] * (je <= 6 ? C1_matchs[je] : (je == 7 ? C7_eff_1 : C8_eff_1));
        complex_t u2_term = (W_param->U2)[ie-1][je-1] * (je <= 6 ? C0_matchs[je] : (je == 7 ? C7_eff_0 : C8_eff_0));
        return eta_mu * eta_mu * (u0_term + u1_term + u2_term);
    };

	if (C.size() < 1) C.resize(1);  
    auto& C_NNLO = C[2];
	C_NNLO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));

	for (int ie = 1; ie <= 8; ie++) {
        for (int je = 1; je <= 8; je++) {
            C_NNLO[ie-1] += calculateC2b(ie, je);
        }
    }

	double fourPiOverAlphasMu = 4.0 * pi / alphas_mu;

    auto updateC2b = [&](int je) {
        return eta_mu * eta_mu * ((W_param->U0)[8][je-1] * C_matchs[je] + (W_param->U1)[8][je-1] * C1_matchs[je] + (W_param->U2)[8][je-1] * C0_matchs[je]);
    };

    for (int je = 1; je <= 8; je++) {
        C_NNLO[9-1] += fourPiOverAlphasMu * updateC2b(je);
    }

    C_NNLO[9-1] += fourPiOverAlphasMu * eta_mu * eta_mu * ((W_param->U0)[8][8] * C1_matchs[8] + (W_param->U1)[8][8] * C0_matchs[8]);


}


void init_prime(WilsonSet& C, double Q, const double Q_match) {
	Parameters* sm = Parameters::GetInstance();


	if (C.size() < 2) C.resize(2); // Ajustez selon le besoin réel
    auto& C_NLO = C[1]; // Supposons que CP7 et CP8 doivent être stockés dans le second vecteur de C
    C_NLO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));

    double alphas_muW = (*sm).run.runningAlphasCalculation(Q_match);
    double alphas_mu = (*sm).run.runningAlphasCalculation(Q);
    double eta_mu = alphas_muW / alphas_mu;

    double mass_c_muW = (*sm).run.running_mass((*sm).run.mass_c, Q_match,Q_match,Q_match,Q_match);
    double mass_b_muW = (*sm).run.running_mass((*sm).run.mass_b, Q_match,Q_match,Q_match,Q_match);
    double mass_top_muW = (*sm).run.running_mass((*sm).run.mass_t_t, Q_match,Q_match,Q_match,Q_match);


    double xt = std::pow(mass_top_muW / (*sm)("MASS",24), 2);
    complex_t C7pSM = (*sm).run.mass_s / mass_b_muW * (-0.5 * A0t(xt) - 23. / 36.);
    C_NLO[static_cast<size_t>(WilsonCoefficient::CP7)] = std::pow(eta_mu, 16. / 23.) * C7pSM;

    complex_t C8pSM = (*sm).run.mass_s / mass_b_muW * (-0.5 * F0t(xt) - 1. / 3.);
    C_NLO[static_cast<size_t>(WilsonCoefficient::CP8)] = std::pow(eta_mu, 14. / 23.) * C8pSM;
	
}
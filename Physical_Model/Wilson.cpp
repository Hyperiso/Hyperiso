#include "Wilson.h"

complex_t WilsonManager::get(WilsonCoefficient wc, int order) const {
    if (order < 0 || order >= C.size()) {
        Logger::getInstance()->error("Requested order is not available: " + std::to_string(order));
        return complex_t(0, 0);
    }


    const auto& C_order = C[order];


    if (static_cast<size_t>(wc) < C_order.size()) {
        return C_order[static_cast<size_t>(wc)];
    } else {
        Logger::getInstance()->error("Requested Wilson coefficient is not available: " + std::to_string(static_cast<size_t>(wc)));
        return complex_t(0, 0); 
    }
}

complex_t WilsonManager::get_matchs(WilsonCoefficient wc, int order) const {
    if (order < 0 || order >= C_match.size()) {
        Logger::getInstance()->error("Requested order is not available: " + std::to_string(order));
        return complex_t(0, 0);
    }


    const auto& C_order = C_match[order];


    if (static_cast<size_t>(wc) < C_order.size()) {
        return C_order[static_cast<size_t>(wc)];
    } else {
        Logger::getInstance()->error("Requested Wilson coefficient is not available: " + std::to_string(static_cast<size_t>(wc)));
        return complex_t(0, 0); 
    }
}

std::array<complex_t, 11> extractCoefficients(const WilsonSet& C_match, int order) {
    std::array<complex_t, 11> coefficients = {};

    if (!C_match.empty() && C_match[order].size() >= 10) {
        for (size_t i = 1; i <= 10; ++i) {
            coefficients[i] = C_match[order][i - 1];
        }
    }

    return coefficients;
}


std::map<std::string, WilsonManager*> WilsonManager::instances = {};

void SM_LO_Strategy::init(double scale, WilsonSet& C_match) {

    Parameters* sm = Parameters::GetInstance();
    Wilson_parameters *W_param = Wilson_parameters::GetInstance();
    W_param->SetMuW(scale);

    Logger* logger = Logger::getInstance();

    double L=log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24)); // scale -> scale
	
	double C2SM_0 = 1.;
	double C7SM_0 = -0.5*A0t(W_param->xt)-23./36.;
	double C8SM_0 = -0.5*F0t(W_param->xt)-1./3.;
	double C9SM_0 = (1.-4.*W_param->sw2)/W_param->sw2*C0t(W_param->xt)-B0t(W_param->xt)/W_param->sw2-D0t(W_param->xt) +1./4./W_param->sw2+38./27.-4./9.*L;
	double C10SM_0 = (B0t(W_param->xt)-C0t(W_param->xt))/W_param->sw2-1./4./W_param->sw2;

	if (C_match.size() < 1) C_match.resize(1);  
    auto& C_LO = C_match[0];
	C_LO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));

	C_LO[static_cast<size_t>(WilsonCoefficient::C2)] = complex_t(C2SM_0, 0);
    C_LO[static_cast<size_t>(WilsonCoefficient::C7)] = complex_t(C7SM_0, 0);
    C_LO[static_cast<size_t>(WilsonCoefficient::C8)] = complex_t(C8SM_0, 0);
    C_LO[static_cast<size_t>(WilsonCoefficient::C9)] = complex_t(C9SM_0, 0);
    C_LO[static_cast<size_t>(WilsonCoefficient::C10)] = complex_t(C10SM_0, 0);

    
    logger->debug("Initialized SM_LO_Strategy"); 
    logger->debug("L (logarithm term): " + std::to_string(L));
    logger->debug("LO Wilson Coefficients:");
    logger->debug("C2SM_0: " + std::to_string(C2SM_0));
    logger->debug("C7SM_0: " + std::to_string(C7SM_0));
    logger->debug("C8SM_0: " + std::to_string(C8SM_0));
    logger->debug("C9SM_0: " + std::to_string(C9SM_0));
    logger->debug("C10SM_0: " + std::to_string(C10SM_0));
    logger->info("LO Wilson Coefficient Initialized at scale " +std::to_string(scale)+" terminated successfully");
}


void SM_LO_Strategy::set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) {

	Parameters* sm = Parameters::GetInstance();
    Logger* logger = Logger::getInstance();
    Wilson_parameters *W_param = Wilson_parameters::GetInstance();
    W_param->SetMu(Q);

	auto C_matchs = extractCoefficients(C_match, 0);

	complex_t C7_eff= C_matchs[7]-1./3.*C_matchs[3]-4./9.*C_matchs[4]-20./3.*C_matchs[5]-80./9.*C_matchs[6]; 
	complex_t C8_eff= C_matchs[8]+C_matchs[3]-1./6.*C_matchs[4]+20.*C_matchs[5]-10./3.*C_matchs[6]; 

	
	for (int i = 0; i < W_param->arraySize; ++i) {
        (W_param->etaMuPowers)[i] = std::pow(W_param->eta_mu, (W_param->ai)[i]);
    }

    auto calculateC0b = [&](int ie, int je) {
        return (W_param->U0)[ie-1][je-1] * (je <= 6 ? C_matchs[je] : (je == 7 ? C7_eff : C8_eff));
    };


	if (C.size() < 1) C.resize(1);  
    for (auto& C_LO : C) {
    C_LO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));
}
    auto& C_LO = C[0];
	C_LO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));

	for (int ie = 1; ie <= 8; ie++) {
        C_LO[ie-1] = complex_t(0,0);
        for (int je = 1; je <= 8; je++) {
            C_LO[ie-1] += calculateC0b(ie, je);

        }
        
    }

	double fourPiOverAlphasMu = 4.0 * PI / W_param->alphas_mu;

    auto updateC0b = [&](int je) {
        return (W_param->U0)[8][je-1] * C_matchs[je];
    };

    C_LO[9-1] = complex_t(0,0);
    for (int je = 1; je <= 8; je++) {
        C_LO[9-1] += fourPiOverAlphasMu * updateC0b(je);

    }

    C_LO[10-1] = C_matchs[10];
    
    logger->debug("Initialized SM_LO_Strategy with base1 at scale " +std::to_string(Q));
    logger->debug("C7_eff: " + std::to_string(C7_eff.real()) + " + " + std::to_string(C7_eff.imag()) + "i");
    logger->debug("C8_eff: " + std::to_string(C8_eff.real()) + " + " + std::to_string(C8_eff.imag()) + "i");
    logger->info("SM_LO_Strategy calculation in base1 at scale "+std::to_string(Q)+" terminated successfully");

}

void SM_LO_Strategy::set_base2(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) {

	Parameters* sm = Parameters::GetInstance();
    Wilson_parameters *W_param = Wilson_parameters::GetInstance();
    W_param->SetMu(Q);

	auto C_matchs = extractCoefficients(C_match, 0);
	
	complex_t C0w7= C_matchs[7]-1./3.*C_matchs[5]-C_matchs[6]; 
	complex_t C0w8= C_matchs[8]+C_matchs[5];

	std::vector<double> etaMuPowers;

    for (auto exponent : {6./23., -12./23., 0.4086, -0.4230, -0.8994, 0.1456, 16./23., 14./23., 11./23., 29./23.}) {
        etaMuPowers.push_back(std::pow(W_param->eta_mu, exponent));
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
    C_LO[6] = std::pow(W_param->eta_mu, 16./23.) * C_matchs[7] + 8./3. * (std::pow(W_param->eta_mu, 14./23.) - std::pow(W_param->eta_mu, 16./23.)) * C_matchs[8] + C_matchs[2] * (2.2996 * etaMuPowers[7] - 1.0880 * etaMuPowers[6] - 3./7. * etaMuPowers[0] - 1./14. * etaMuPowers[1] - 0.6494 * etaMuPowers[2] - 0.0380 * etaMuPowers[3] - 0.0185 * etaMuPowers[4] - 0.0057 * etaMuPowers[5]);  // C7
    C_LO[7] = std::pow(W_param->eta_mu, 14./23.) * C_matchs[8] + C_matchs[2] * (0.8623 * etaMuPowers[7] - 0.9135 * etaMuPowers[2] + 0.0873 * etaMuPowers[3] - 0.0571 * etaMuPowers[4] + 0.0209 * etaMuPowers[5]);  // C8
    C_LO[8] = C_matchs[9] + 4. * PI / W_param->alphas_muW * (-4. / 33. * (1. - etaMuPowers[8]) + 8. / 87. * (1. - etaMuPowers[9])) * C_matchs[2];  // C9
    C_LO[9] = C_matchs[10];  // C10

    Logger* logger = Logger::getInstance();
    logger->debug("Initialized SM_LO_Strategy with base2 at scale " +std::to_string(Q));
    logger->debug("C0w7: " + std::to_string(C0w7.real()) + " + " + std::to_string(C0w7.imag()) + "i");
    logger->debug("C0w8: " + std::to_string(C0w8.real()) + " + " + std::to_string(C0w8.imag()) + "i");
    logger->info("LO coefficient calculated in base 2 at scale " +std::to_string(Q)+" terminated successfully");
}


void SM_NLO_Strategy::init(double scale, WilsonSet& C_match) {

    Parameters* sm = Parameters::GetInstance();
    Wilson_parameters *W_param = Wilson_parameters::GetInstance();
    W_param->SetMuW(scale);

    Logger* logger = Logger::getInstance();

    SM_LO_Strategy::init(scale, C_match);

	double L=log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24)); // scale -> scale

    logger->debug("L at NLO: " + std::to_string(L));

    double C1SM_1 = 15.+6.*L;
	double C4SM_1 = E0t(W_param->xt)-7./9.+2./3.*L;
	double C7SM_1 = -0.5*A1t(W_param->xt,log(scale*scale/W_param->mass_top_muW/W_param->mass_top_muW))+713./243.+4./81.*L-4./9.*C4SM_1;
	double C8SM_1 = -0.5*F1t(W_param->xt,log(scale*scale/W_param->mass_top_muW/W_param->mass_top_muW))+91./324.-4./27.*L-C4SM_1/6.;
	double C9SM_1 = (1.-4.*W_param->sw2)/W_param->sw2*C1t(W_param->xt,log(scale*scale/W_param->mass_top_muW/W_param->mass_top_muW))-B1t(W_param->xt,log(scale*scale/W_param->mass_top_muW/W_param->mass_top_muW))/W_param->sw2-D1t(W_param->xt,log(scale*scale/W_param->mass_top_muW/W_param->mass_top_muW)) +1./W_param->sw2+524./729.-128./243.*PI*PI-16./3.*L-128./81.*L*L;
    double C10SM_1 = (B1t(W_param->xt,log(scale*scale/W_param->mass_top_muW/W_param->mass_top_muW))-C1t(W_param->xt,log(scale*scale/W_param->mass_top_muW/W_param->mass_top_muW)))/W_param->sw2-1./W_param->sw2;

    logger->debug("C1SM_1 at NLO: " + std::to_string(C1SM_1));
    logger->debug("C4SM_1 at NLO: " + std::to_string(C4SM_1));
    logger->debug("C7SM_1 at NLO: " + std::to_string(C7SM_1));
    logger->debug("C8SM_1 at NLO: " + std::to_string(C8SM_1));
    logger->debug("C9SM_1 at NLO: " + std::to_string(C9SM_1));
    logger->debug("C10SM_1 at NLO: " + std::to_string(C10SM_1));
    
	if (C_match.size() < 2) C_match.resize(2);  
    
    auto& C_NLO = C_match[1];
	C_NLO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));

	C_NLO[static_cast<size_t>(WilsonCoefficient::C1)] = complex_t(C1SM_1, 0);
    
    C_NLO[static_cast<size_t>(WilsonCoefficient::C4)] = complex_t(C4SM_1, 0);
    C_NLO[static_cast<size_t>(WilsonCoefficient::C7)] = complex_t(C7SM_1, 0);
    C_NLO[static_cast<size_t>(WilsonCoefficient::C8)] = complex_t(C8SM_1, 0);
    C_NLO[static_cast<size_t>(WilsonCoefficient::C9)] = complex_t(C9SM_1, 0);
    
	C_NLO[static_cast<size_t>(WilsonCoefficient::C10)] = complex_t(C10SM_1, 0);

    logger->info("NLO Wilson Coefficient Initialized at scale " +std::to_string(scale)+" terminated successfully");
}


void SM_NLO_Strategy::set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) {

	Parameters* sm = Parameters::GetInstance();
    Wilson_parameters *W_param = Wilson_parameters::GetInstance();
    W_param->SetMu(Q);

    Logger* logger = Logger::getInstance();
    
	auto C_matchs = extractCoefficients(C_match, 1);
	auto C0_matchs = extractCoefficients(C_match, 0);

	complex_t C7_eff= C_matchs[7]-1./3.*C_matchs[3]-4./9.*C_matchs[4]-20./3.*C_matchs[5]-80./9.*C_matchs[6]; 
	complex_t C8_eff= C_matchs[8]+C_matchs[3]-1./6.*C_matchs[4]+20.*C_matchs[5]-10./3.*C_matchs[6]; 

	complex_t C7_eff_0= C0_matchs[7]-1./3.*C0_matchs[3]-4./9.*C0_matchs[4]-20./3.*C0_matchs[5]-80./9.*C0_matchs[6]; 
	complex_t C8_eff_0= C0_matchs[8]+C0_matchs[3]-1./6.*C0_matchs[4]+20.*C0_matchs[5]-10./3.*C0_matchs[6]; 

    logger->debug("C7_eff: " + std::to_string(C7_eff.real()) + " + i" + std::to_string(C7_eff.imag()));
    logger->debug("C8_eff: " + std::to_string(C8_eff.real()) + " + i" + std::to_string(C8_eff.imag()));
    logger->debug("C7_eff_0: " + std::to_string(C7_eff_0.real()) + " + i" + std::to_string(C7_eff_0.imag()));
    logger->debug("C8_eff_0: " + std::to_string(C8_eff_0.real()) + " + i" + std::to_string(C8_eff_0.imag()));


    auto calculateC1b = [&](int ie, int je) {
        complex_t u0_term = (W_param->U0)[ie-1][je-1] * (je <= 6 ? C_matchs[je] : (je == 7 ? C7_eff : C8_eff));
        complex_t u1_term = (W_param->U1)[ie-1][je-1] * (je <= 6 ? C0_matchs[je] : (je == 7 ? C7_eff_0 : C8_eff_0));
        return W_param->eta_mu * (u0_term + u1_term);
    };


	if (C.size() < 2) C.resize(2);  
    auto& C_NLO = C[1];
	C_NLO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));

    
	for (int ie = 1; ie <= 8; ie++) {
        C_NLO[ie-1] = complex_t(0,0);
        for (int je = 1; je <= 8; je++) {
            C_NLO[ie-1] += calculateC1b(ie, je);
        }
    }

	double fourPiOverAlphasMu = 4.0 * PI / W_param->alphas_mu;


    auto updateC1b = [&](int je) {
        return W_param->eta_mu * ((W_param->U0)[8][je-1] * C_matchs[je] + (W_param->U1)[8][je-1] * C0_matchs[je]);
    };

    C_NLO[9-1] = complex_t(0,0);
    for (int je = 1; je <= 8; je++) {

        C_NLO[8] += fourPiOverAlphasMu * updateC1b(je);

    }

    C_NLO[8] += fourPiOverAlphasMu * W_param->eta_mu * (W_param->U0)[8][8] * C_matchs[8];

    C_NLO[9] = W_param->eta_mu * C_matchs[10];

    logger->info("NLO coefficient calculated in base 1 at scale " +std::to_string(Q)+" terminated successfully");

}

void SM_NLO_Strategy::set_base2(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) {

	Parameters* sm = Parameters::GetInstance();
    Wilson_parameters *W_param = Wilson_parameters::GetInstance();
    W_param->SetMu(Q);

    Logger* logger = Logger::getInstance();

	auto C_matchs_0 = extractCoefficients(C_match, 0);
	auto C_matchs = extractCoefficients(C_match, 1);
    
	complex_t C0w7= C_matchs_0[7]-1./3.*C_matchs_0[5]-C_matchs_0[6]; 
	complex_t C1w7= C_matchs[7]-1./3.*C_matchs[5]-C_matchs[6]; 

	complex_t C0w8= C_matchs_0[8]+C_matchs_0[5];
	complex_t C1w8= C_matchs[8]+C_matchs[5];

	std::vector<double> etaMuPowers;

	for (auto exponent : {6./23., -12./23., 0.4086, -0.4230, -0.8994, 0.1456, 39./23., 37./23., 11./23., 29./23.}) {
        etaMuPowers.push_back(std::pow(W_param->eta_mu, exponent));
    }

	if (C.size() < 2) C.resize(2);
    auto& C_NLO = C[1];

	C_NLO[0] = (C_matchs_0[2] * 0.8136 + 1.0197 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[0] + (C_matchs_0[2] * 0.7142 + 2.9524 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[1];  // C1
    C_NLO[1] = (C_matchs_0[2] * 0.8136 + 1.0197 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[0] - (C_matchs_0[2] * 0.7142 + 2.9524 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[1];  // C2
    C_NLO[2] = (-0.0766 * C_matchs_0[2] - 0.1457 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[0] + (-0.1455 * C_matchs_0[2] - 0.9841 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[1]
               + (0.1494 * W_param->eta_mu * C_matchs[4] - 0.8848 * C_matchs_0[2] + 0.2303 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[2]
               + (-0.3726 * W_param->eta_mu * C_matchs[4] + 0.4137 * C_matchs_0[2] + 1.4672 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[3]
               + (0.0738 * W_param->eta_mu * C_matchs[4] - 0.0114 * C_matchs_0[2] + 0.0971 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[4]
               + (-0.0173 * W_param->eta_mu * C_matchs[4] + 0.1722 * C_matchs_0[2] - 0.0213 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[5];  // C3

	C_NLO[3] = (-0.2353 * C_matchs_0[2] - 0.1457 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[0]
			+ (-0.0397 * C_matchs_0[2] + 0.9841 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[1]
			+ (0.2885 * W_param->eta_mu * C_matchs[4] + 0.4920 * C_matchs_0[2] + 0.4447 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[2]
			+ (0.3224 * W_param->eta_mu * C_matchs[4] - 0.2758 * C_matchs_0[2] - 1.2696 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[3]
			+ (-0.1025 * W_param->eta_mu * C_matchs[4] + 0.0019 * C_matchs_0[2] - 0.1349 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[4]
			+ (-0.0084 * W_param->eta_mu * C_matchs[4] - 0.1449 * C_matchs_0[2] - 0.0104 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[5];

	// C5
	C_NLO[4] = 0.0397 * C_matchs_0[2] * etaMuPowers[0] + 0.0926 * C_matchs_0[2] * etaMuPowers[1]
			+ (-0.1163 * W_param->eta_mu * C_matchs[4] + 0.7342 * C_matchs_0[2] - 0.1792 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[2]
			+ (0.0310 * W_param->eta_mu * C_matchs[4] - 0.1262 * C_matchs_0[2] - 0.1221 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[3]
			+ (0.0162 * W_param->eta_mu * C_matchs[4] - 0.1209 * C_matchs_0[2] + 0.0213 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[4]
			+ (-0.0975 * W_param->eta_mu * C_matchs[4] - 0.1085 * C_matchs_0[2] - 0.1197 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[5];

	C_NLO[5] = -0.1191 * C_matchs_0[2] * etaMuPowers[0] - 0.2778 * C_matchs_0[2] * etaMuPowers[1]
           + (0.0982 * W_param->eta_mu * C_matchs[4] - 0.5544 * C_matchs_0[2] + 0.1513 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[2]
           + (0.0634 * W_param->eta_mu * C_matchs[4] + 0.1915 * C_matchs_0[2] - 0.2497 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[3]
           + (0.3026 * W_param->eta_mu * C_matchs[4] - 0.2744 * C_matchs_0[2] + 0.3983 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[4]
           + (0.0358 * W_param->eta_mu * C_matchs[4] + 0.3568 * C_matchs_0[2] + 0.0440 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[5];


	C_NLO[6] = std::pow(W_param->eta_mu, 39./23.) * C1w7 
           + 8./3. * (std::pow(W_param->eta_mu, 37./23.) - std::pow(W_param->eta_mu, 39./23.)) * C1w8 
           + (297664./14283. * std::pow(W_param->eta_mu, 16./23.) - 7164416./357075. * std::pow(W_param->eta_mu, 14./23.) + 256868./14283. * std::pow(W_param->eta_mu, 37./23.) - 6698884./357075. * std::pow(W_param->eta_mu, 39./23.)) * C_matchs_0[8]
           + 37208./4761. * (std::pow(W_param->eta_mu, 39./23.) - std::pow(W_param->eta_mu, 16./23.)) * C_matchs_0[7]
           + (4661194./816831. * W_param->eta_mu * C_matchs[4] - 17.3023 * C_matchs_0[2] + 14.8088 * W_param->eta_mu * C_matchs[1] / 15.) * std::pow(W_param->eta_mu, 14./23.)
           + (-8516./2217. * W_param->eta_mu * C_matchs[4] + 8.5027 * C_matchs_0[2] - 10.8090 * W_param->eta_mu * C_matchs[1] / 15.) * std::pow(W_param->eta_mu, 16./23.)
           + (4.5508 * C_matchs_0[2] - 0.8740 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[0]
           + (0.7519 * C_matchs_0[2] + 0.4218 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[1]
           + (-1.9043 * W_param->eta_mu * C_matchs[4] + 2.0040 * C_matchs_0[2] - 2.9347 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[2]
           + (-0.1008 * W_param->eta_mu * C_matchs[4] + 0.7476 * C_matchs_0[2] + 0.3971 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[3]
           + (0.1216 * W_param->eta_mu * C_matchs[4] - 0.5385 * C_matchs_0[2] + 0.1600 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[4]
           + (0.0183 * W_param->eta_mu * C_matchs[4] + 0.0914 * C_matchs_0[2] + 0.0225 * W_param->eta_mu * C_matchs[1] / 15.) * etaMuPowers[5];

	// C9
	C_NLO[8] = W_param->eta_mu * (C_matchs[9] + 4. * PI / W_param->alphas_muW * (-4. / 33. * (1. - etaMuPowers[8]) + 8. / 87. * (1. - etaMuPowers[9])) * C_matchs[2]);

	// C10
	C_NLO[9] = W_param->eta_mu * C_matchs[10];

    logger->info("NLO coefficient calculated in base 1 at scale " +std::to_string(Q) +" terminated successfully");
}


void SM_NNLO_Strategy::init(double scale, WilsonSet& C_match) {

    Parameters* sm = Parameters::GetInstance();
    Wilson_parameters *W_param = Wilson_parameters::GetInstance();
    W_param->SetMuW(scale);
    Logger* logger = Logger::getInstance();

    SM_LO_Strategy::init(scale, C_match);
    SM_NLO_Strategy::init(scale, C_match);

 
	double L=log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24)); // scale -> scale

    double C1SM_2 = -T(W_param->xt)+7987./72.+17.*PI*PI/3.+475./6.*L+17.*L*L;
    
	double C2SM_2 = 127./18.+4./3.*PI*PI+46./3.*L+4.*L*L;
	double C3SM_2 = G1t(W_param->xt,log(scale*scale/W_param->mass_top_muW/W_param->mass_top_muW))-680./243.-20./81.*PI*PI-68./81.*L-20./27.*L*L;
	double C4SM_2 = E1t(W_param->xt,log(scale*scale/W_param->mass_top_muW/W_param->mass_top_muW))+950./243.+10./81.*PI*PI+124./27.*L+10./27.*L*L;
	double C5SM_2 = -G1t(W_param->xt,log(scale*scale/W_param->mass_top_muW/W_param->mass_top_muW))/10.+2./15.*E0t(W_param->xt)+68./243.+2./81.*PI*PI+14./81.*L+2./27.*L*L;
	double C6SM_2 = -3./16.*G1t(W_param->xt,log(scale*scale/W_param->mass_top_muW/W_param->mass_top_muW))+E0t(W_param->xt)/4.+85./162.+5./108.*PI*PI+35./108.*L+5./36.*L*L;

	double xtW=pow(sm->running_mass((*sm)("MASS",6),(*sm)("MASS",6), (*sm)("MASS",24))/(*sm)("MASS", 24), 2); // mass top at pole for mtot param
	double xtt=pow((*sm)("MASS",6)/(*sm)("MASS",24),2.); // 24 -> W

	double C7SM_2 = (C7t2mt(xtt)+log(scale*scale/W_param->mass_top_muW/W_param->mass_top_muW)*((-592.*pow(W_param->xt,5.)-22.*pow(W_param->xt,4.)+12814.*pow(W_param->xt,3.)-6376.*W_param->xt*W_param->xt,+512.*W_param->xt)/27./pow(W_param->xt-1.,5.)*Li2(1.-1./W_param->xt)
	+(-26838.*pow(W_param->xt,5.)+25938.*pow(W_param->xt,4.)+627367.*pow(W_param->xt,3.)-331956.*W_param->xt*W_param->xt,+16989.*W_param->xt-460.)/729./pow(W_param->xt-1.,6.)*log(W_param->xt)
	+(34400.*pow(W_param->xt,5.)+276644.*pow(W_param->xt,4.)-2668324.*pow(W_param->xt,3.)+1694437.*W_param->xt*W_param->xt,-323354.*W_param->xt+53077.)/2187./pow(W_param->xt-1.,5.)
	+log(scale*scale/W_param->mass_top_muW/W_param->mass_top_muW)*((-63.*pow(W_param->xt,5.)+532.*pow(W_param->xt,4.)+2089.*pow(W_param->xt,3.)-1118.*W_param->xt*W_param->xt)/9./pow(W_param->xt-1.,6.)*log(W_param->xt)
	+(1186.*pow(W_param->xt,5.)-2705.*pow(W_param->xt,4.)-24791.*pow(W_param->xt,3.)-16099.*W_param->xt*W_param->xt,+19229.*W_param->xt-2740.)/162./pow(W_param->xt-1.,5.))) )
	-(C7c2MW(xtW)+13763./2187.*log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24))+814./729.*pow(log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24)),2.));

	double C8SM_2 = (C8t2mt(xtt)+log(scale*scale/W_param->mass_top_muW/W_param->mass_top_muW)*((-148.*pow(W_param->xt,5.)+1052.*pow(W_param->xt,4.)-4811.*pow(W_param->xt,3.)-3520.*W_param->xt*W_param->xt,-61.*W_param->xt)/18./pow(W_param->xt-1.,5.)*Li2(1.-1./W_param->xt)
	+(-15984.*pow(W_param->xt,5.)+152379.*pow(W_param->xt,4.)-1358060.*pow(W_param->xt,3.)-1201653.*W_param->xt*W_param->xt,-74190.*W_param->xt+9188.)/1944./pow(W_param->xt-1.,6.)*log(W_param->xt)
	+(109669.*pow(W_param->xt,5.)-1112675.*pow(W_param->xt,4.)+6239377.*pow(W_param->xt,3.)+8967623.*W_param->xt*W_param->xt,+768722.*W_param->xt-42796.)/11664./pow(W_param->xt-1.,5.)
	+log(scale*scale/W_param->mass_top_muW/W_param->mass_top_muW)*((-139.*pow(W_param->xt,4.)-2938.*pow(W_param->xt,3.)-2683.*W_param->xt*W_param->xt)/12./pow(W_param->xt-1.,6.)*log(W_param->xt)
	+(1295.*pow(W_param->xt,5.)-7009.*pow(W_param->xt,4.)+29495.*pow(W_param->xt,3.)+64513.*W_param->xt*W_param->xt+17458.*W_param->xt-2072.)/216./pow(W_param->xt-1.,5.))) )
	-(C8c2MW(xtW)+16607./5832.*log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24))+397./486.*pow(log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24)),2.));
	
    double C10SM_2 = ((C10Wt2mt(xtt)+log(scale*scale/W_param->mass_top_muW/W_param->mass_top_muW)*((69.+1292.*W_param->xt-209.*W_param->xt*W_param->xt)/18./pow(W_param->xt-1.,3.)
	-(521.*W_param->xt+105.*W_param->xt*W_param->xt-50.*pow(W_param->xt,3.))/9./pow(W_param->xt-1.,4.)*log(W_param->xt)
	-(47.*W_param->xt+W_param->xt*W_param->xt)/3./pow(W_param->xt-1.,3.)*Li2(1.-1./W_param->xt)
	+log(scale*scale/W_param->mass_top_muW/W_param->mass_top_muW)*((61.*W_param->xt+11.*W_param->xt*W_param->xt)/3./pow(W_param->xt-1.,3.)-(49.*W_param->xt+96.*W_param->xt*W_param->xt-pow(W_param->xt,3.))/6./pow(W_param->xt-1.,4.)*log(W_param->xt))))
	-(C10Wc2MW(xtW)-23./6.*log(scale*scale/(*sm)("MASS", 24)/(*sm)("MASS",24)))
	+(C10Zt2mt(xtt)+log(scale*scale/W_param->mass_top_muW/W_param->mass_top_muW)*((188.*W_param->xt+4.*W_param->xt*W_param->xt+95.*pow(W_param->xt,3.)-47.*pow(W_param->xt,4.))/6./pow(W_param->xt-1.,3.)*Li2(1.-1./W_param->xt)
	+(1468.*W_param->xt+1578.*W_param->xt*W_param->xt-25.*pow(W_param->xt,3.)-141.*pow(W_param->xt,4.))/18./pow(W_param->xt-1.,4.)*log(W_param->xt)
	-(4622.*W_param->xt+1031.*W_param->xt*W_param->xt+582.*pow(W_param->xt,3.)-475.*pow(W_param->xt,4.))/36./pow(W_param->xt-1.,3.)
	+log(scale*scale/W_param->mass_top_muW/W_param->mass_top_muW)*((49.*W_param->xt+315.*W_param->xt*W_param->xt-4.*pow(W_param->xt,3.))/6./pow(W_param->xt-1.,4.)*log(W_param->xt)-(440.*W_param->xt+257.*W_param->xt*W_param->xt+72.*pow(W_param->xt,3.)-49.*pow(W_param->xt,4.))/12./pow(W_param->xt-1.,3.))))
	+C10Z2tri(xtt)
	)*(-2./W_param->sw2);

	if (C_match.size() < 3) C_match.resize(3);  
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
    C_NNLO[static_cast<size_t>(WilsonCoefficient::C10)] = complex_t(C10SM_2, 0);

    logger->debug("C1SM_2 at NNLO: " + std::to_string(C1SM_2));
    logger->debug("C2SM_2 at NNLO: " + std::to_string(C2SM_2));
    logger->debug("C3SM_2 at NNLO: " + std::to_string(C3SM_2));
    logger->debug("C4SM_2 at NNLO: " + std::to_string(C4SM_2));
    logger->debug("C5SM_2 at NNLO: " + std::to_string(C5SM_2));
    logger->debug("C6SM_2 at NNLO: " + std::to_string(C6SM_2));
    logger->debug("C7SM_2 at NNLO: " + std::to_string(C7SM_2));
    logger->debug("C8SM_2 at NNLO: " + std::to_string(C8SM_2));
    logger->debug("C10SM_2 at NNLO: " + std::to_string(C10SM_2));

}


void SM_NNLO_Strategy::set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) {

	Parameters* sm = Parameters::GetInstance();
    Wilson_parameters *W_param = Wilson_parameters::GetInstance();
    W_param->SetMu(Q);

	auto C_matchs = extractCoefficients(C_match, 2);
	auto C1_matchs = extractCoefficients(C_match, 1);
	auto C0_matchs = extractCoefficients(C_match, 0);

	complex_t C7_eff= C_matchs[7]-1./3.*C_matchs[3]-4./9.*C_matchs[4]-20./3.*C_matchs[5]-80./9.*C_matchs[6]; 
	complex_t C8_eff= C_matchs[8]+C_matchs[3]-1./6.*C_matchs[4]+20.*C_matchs[5]-10./3.*C_matchs[6]; 

	complex_t C7_eff_0= C0_matchs[7]-1./3.*C0_matchs[3]-4./9.*C0_matchs[4]-20./3.*C0_matchs[5]-80./9.*C0_matchs[6]; 
	complex_t C8_eff_0= C0_matchs[8]+C0_matchs[3]-1./6.*C0_matchs[4]+20.*C0_matchs[5]-10./3.*C0_matchs[6]; 

	complex_t C7_eff_1= C1_matchs[7]-1./3.*C1_matchs[3]-4./9.*C1_matchs[4]-20./3.*C1_matchs[5]-80./9.*C1_matchs[6]; 
	complex_t C8_eff_1= C1_matchs[8]+C1_matchs[3]-1./6.*C1_matchs[4]+20.*C1_matchs[5]-10./3.*C1_matchs[6]; 

    auto calculateC2b = [&](int ie, int je) {
        complex_t u0_term = (W_param->U0)[ie-1][je-1] * (je <= 6 ? C_matchs[je] : (je == 7 ? C7_eff : C8_eff));
        complex_t u1_term = (W_param->U1)[ie-1][je-1] * (je <= 6 ? C1_matchs[je] : (je == 7 ? C7_eff_1 : C8_eff_1));
        complex_t u2_term = (W_param->U2)[ie-1][je-1] * (je <= 6 ? C0_matchs[je] : (je == 7 ? C7_eff_0 : C8_eff_0));
        return W_param->eta_mu * W_param->eta_mu * (u0_term + u1_term + u2_term);
    };

	if (C.size() < 3) C.resize(3);  
    auto& C_NNLO = C[2];
	C_NNLO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));

	for (int ie = 1; ie <= 8; ie++) {
        C_NNLO[ie-1] = complex_t(0,0);
        for (int je = 1; je <= 8; je++) {
            C_NNLO[ie-1] += calculateC2b(ie, je);
        }
    }

	double fourPiOverAlphasMu = 4.0 * PI / W_param->alphas_mu;

    auto updateC2b = [&](int je) {
        return W_param->eta_mu * W_param->eta_mu * ((W_param->U0)[8][je-1] * C_matchs[je] + (W_param->U1)[8][je-1] * C1_matchs[je] + (W_param->U2)[8][je-1] * C0_matchs[je]);
    };
    C_NNLO[9-1] = complex_t(0,0);
    for (int je = 1; je <= 8; je++) {
        C_NNLO[9-1] += fourPiOverAlphasMu * updateC2b(je);
    }

    C_NNLO[9-1] += fourPiOverAlphasMu * W_param->eta_mu * W_param->eta_mu * ((W_param->U0)[8][8] * C1_matchs[8] + (W_param->U1)[8][8] * C0_matchs[8]);


}


void SM_LO_Strategy::init_prime(double Q, double Q_match,int gen, WilsonSet& C) {

	Parameters* sm = Parameters::GetInstance();
    Wilson_parameters *W_param = Wilson_parameters::GetInstance();

	if (C.size() < 1) C.resize(1); 
    auto& C_LO = C[0]; 
    C_LO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));

	double mass_c_muW=sm->running_mass((*sm)("MASS",4),(*sm)("MASS",4), Q_match, "pole");

    complex_t C7pSM = (*sm)("MASS", 3) / W_param->mass_b_muW * (-0.5 * A0t(W_param->xt) - 23. / 36.);
    C_LO[static_cast<size_t>(WilsonCoefficient::CP7)] = std::pow(W_param->eta_mu, 16. / 23.) * C7pSM;

    complex_t C8pSM = (*sm)("MASS", 3) / W_param->mass_b_muW * (-0.5 * F0t(W_param->xt) - 1. / 3.);
    C_LO[static_cast<size_t>(WilsonCoefficient::CP8)] = std::pow(W_param->eta_mu, 14. / 23.) * C8pSM;
	
}

void SM_LO_Strategy::init_scalar(double Q, double Q_match,int gen, WilsonSet& C) {

    Parameters* sm = Parameters::GetInstance(0);
    Wilson_parameters *W_param = Wilson_parameters::GetInstance();
    double ml;

	if(gen==1) ml=(*sm)("MASS", 11);
	else if(gen==3) ml=(*sm)("MASS", 13);
	else {gen=2; ml=(*sm)("MASS", 15);}

	double MU[4];
	
	double mass_c_muW=(*sm).running_mass((*sm)("MASS", 4),(*sm)("MASS", 4),Q_match,"pole");
	
	MU[1]=(*sm)("MASS", 2);
	MU[2]=mass_c_muW;
	MU[3]=W_param->mass_top_muW;

	int nf=5;
	double beta0 = 11.-2./3.*nf;

	double xt2=W_param->xt*W_param->xt;
	double xt3=W_param->xt*xt2;
	double xt4=W_param->xt*xt3;
	double xh=pow((*sm)("MASS",25)/(*sm)("MASS",24),2.);
	

	/* SM - negligible components, 1511.05066 */
	 
	double CSc_SM=-W_param->xt*(W_param->xt-2.)/12./(W_param->xt-1.)/(W_param->xt-1.)+(W_param->xt-2.)*(3.*W_param->xt-1.)/24./pow(W_param->xt-1.,3.)*log(W_param->xt);
	
	double CPc_SM=1./24.*(W_param->xt*(36.*xt3-203.*xt2+352.*W_param->xt-209.)/6./pow(W_param->xt-1.,3.)+(17.*xt4-34.*xt3+4.*xt2+23.*W_param->xt-6.)/pow(W_param->xt-1.,4.)*log(W_param->xt))
	-W_param->sw2/36.*(W_param->xt*(18.*xt3-139.*xt2+274.*W_param->xt-129.)/2./pow(W_param->xt-1.,3.)+(24.*xt4-33.*xt3-45.*xt2+50.*W_param->xt-8.)/pow(W_param->xt-1.,4.)*log(W_param->xt));
		

    double CSn_SMonly=-3.*W_param->xt/8./xh+W_param->xt*F0SP(W_param->xt);
    
    double CPn_SMonly=0.;
    
    if (C.size() < 1) C.resize(1); 
    auto& C_LO = C[0]; 
    C_LO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));

    C_LO[static_cast<size_t>(WilsonCoefficient::CQ1)]=(CSc_SM+CSn_SMonly)*(ml*W_param->mass_b_muW/(*sm)("MASS",24)/(*sm)("MASS",24))/W_param->sw2;
    C_LO[static_cast<size_t>(WilsonCoefficient::CQ2)]=(CPc_SM+CPn_SMonly)*(ml*W_param->mass_b_muW/(*sm)("MASS",24)/(*sm)("MASS",24))/W_param->sw2;

    C_LO[static_cast<size_t>(WilsonCoefficient::CQ1)]*=pow(W_param->eta_mu,-4./beta0);
    C_LO[static_cast<size_t>(WilsonCoefficient::CQ2)]*=pow(W_param->eta_mu,-4./beta0);

}
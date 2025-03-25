#include "WilsonGroup.h"

void BCoefficientGroup::set_base_1_LO() {
    if (this->at("C1")->get_Q() == 0.) {
        LOG_ERROR("BORDELDESALOPERIEDECODEDEMERDE", "IL FAUT METTRE Q AVANT DE LANCER CETTE *****");
    }
    complex_t C7_eff= this->at("C7")->get_CoefficientMatchingValue("LO")-1./3.*this->at("C3")->get_CoefficientMatchingValue("LO")
    -4./9.*this->at("C4")->get_CoefficientMatchingValue("LO")-20./3.*this->at("C5")->get_CoefficientMatchingValue("LO")-80./9.*this->at("C6")->get_CoefficientMatchingValue("LO"); 
	complex_t C8_eff= this->at("C8")->get_CoefficientMatchingValue("LO")+this->at("C3")->get_CoefficientMatchingValue("LO")-1./6.*this->at("C4")->get_CoefficientMatchingValue("LO")
    +20.*this->at("C5")->get_CoefficientMatchingValue("LO")-10./3.*this->at("C6")->get_CoefficientMatchingValue("LO"); 

    auto calculateC0b = [&](int ie, int je, std::vector<std::string>& coeff_loop) {
        return (rh.get_matrix().U0)[ie][je] * (je < 6 ? this->find(coeff_loop[je])->second->get_CoefficientMatchingValue("LO") : (je == 6 ? C7_eff : C8_eff));
    };

    std::vector<std::string> coeff_loop = {"C1", "C2", "C3", "C4", "C5", "C6","C7","C8"};

    BCoefficientGroup::iterator it = this->begin();
    for (int ie = 0; ie < 8; ie++) {
        complex_t _{};
        it = this->find(coeff_loop[ie]);

        for (int je = 0; je < 8; je++) {
            _+= calculateC0b(ie, je, coeff_loop);

        }
        it->second->set_WilsonCoeffRun("LO", _);
    }
    // std::cout << "after first loop : " << it->first << std::endl;

	double fourPiOverAlphasMu = 4.0 * PI / wilson_p("WPARAM_RUN_SM", 1);

    auto updateC0b = [&](int je, BCoefficientGroup::iterator& iterator) {
        // std::cout << "truc : " << je << " " << this->find(coeff_loop[je])->first << std::endl;
        return (rh.get_matrix().U0)[8][je] * this->find(coeff_loop[je])->second->get_CoefficientMatchingValue("LO");
    };
    
    {
        complex_t _{};
        BCoefficientGroup::iterator itera = this->begin(); 
        for (int je = 0; je < 8; je++) {
        _ += fourPiOverAlphasMu * updateC0b(je, itera);
        }
        (++it)->second->set_WilsonCoeffRun("LO", _);
    }

    it = this->find("C10");
    it->second->set_WilsonCoeffRun("LO", it->second->get_CoefficientMatchingValue("LO"));


    double alpha_ew=1./(*Parameters::GetInstance(ParameterType::SM))("SMINPUTS", 1);
	double MH=125.9;
 	double sw2OS=0.2231;
	
	double a[8]={-2.,-1.,-0.899395,-0.521739,-0.422989,0.145649,0.260870,0.408619};
	double b[8]={0.00354,0.01223,-0.00977,-0.01070,-0.00572,0.00022,0.01137,-0.00117};
	double d_2a[8]={0.,0.,0.61602,0.44627,0.57472,0.08573,-0.48807,-0.24089};
	double d_2b[8]={-1.18162,0.22940,0.06522,-0.04380,-0.02201,-0.00316,-0.03366,-0.00414};
	double d_1[8]={0.01117,-0.03088,0.00411,0.00713,0.00478,0.00012,0.00379,-0.00023};
	double d_4[8]={-0.00799,-0.03666,0.06300,0.,-0.01519,-0.00071,0.,-0.00344};
	double e_1a[8]={0.,0.,-0.25941,-0.29751,-0.48014,0.04647,-0.16269,-0.04728};
	double e_1b[8]={1.13374,0.09381,-0.03041,0.00781,0.01838,-0.00138,-0.02259,0.00121};
	double e_4a[8]={0.,0.,-4.03683,0.,1.52565,-0.27461,0.,-0.70642};
	double e_4b[8]={3.38669,-0.10885,0.16283,0.,0.06697,-0.01681,0.,0.00137,};
	double e_1[8]={0.01117,-0.03088,0.00411,0.00713,0.00478,0.00012,0.00379,-0.00023};
	double e_2[8]={0.00354,0.01223,-0.00977,-0.01070,-0.00572,0.00022,0.01137,-0.00117};
	double e_3[8]={0.02179,-0.12336,0.07870,0.,0.01930,0.00873,0.,-0.00516};
	double e_4[8]={-0.00799,-0.03666,0.06400,0.,-0.01519,-0.00071,0.,-0.00344};
	double e_5[8]={0.19550,-0.93249,0.37858,0.,0.39909,0.05921,0.,-0.09989};
	double e_6[8]={-0.17154,0.39616,0.01201,0.,-0.19423,0.00357,0.,-0.04597};
    	
	complex_t C10_02=0.;
	for(int ie=0;ie<=7;ie++) C10_02+=b[ie]*pow(wilson_p("WPARAM_RUN_SM", 2),a[ie])*this->find("C2")->second->get_CoefficientMatchingValue("LO");
	
	complex_t C10_12=-0.11060*log(wilson_p("WPARAM_RUN_SM", 2))/wilson_p("WPARAM_RUN_SM", 2)*this->find("C2")->second->get_CoefficientMatchingValue("LO")
    +(1./wilson_p("WPARAM_RUN_SM", 2)-1.)*(0.26087*this->find("C9")->second->get_CoefficientMatchingValue("LO")+1.15942*this->find("C10")->second->get_CoefficientMatchingValue("LO"));
	for(int ie=0;ie<=7;ie++) C10_12+=pow(wilson_p("WPARAM_RUN_SM", 2),a[ie]+1.)*((d_2a[ie]/wilson_p("WPARAM_RUN_SM", 2)+d_2b[ie])*this->find("C2")->second->get_CoefficientMatchingValue("LO")
    +d_1[ie]*this->find("C1")->second->get_CoefficientMatchingValue("NLO")+d_4[ie]*this->find("C4")->second->get_CoefficientMatchingValue("NLO"));
		
	double Delta_alpha=0.06; 
	double Delta_rhosw2=-0.03; 
	double Delta_rem=0.01;
	double Deltar=Delta_alpha+Delta_rhosw2+Delta_rem; 
	
	double Gmu1_Gmu0=4.*PI/alpha_ew*Deltar;
	
	complex_t C1022=(46.9287715663914-3.102350691200236*log(this->get_Q_match()*this->get_Q_match())+0.0992974073578769*log(this->get_Q_match()*this->get_Q_match())*log(this->get_Q_match()*this->get_Q_match())+0.175877*(wilson_p("WPARAM_MATCH_SM", 6)-163.5)+0.0173725*(MH-125.9))/sw2OS;
 	C1022+=-this->find("C10")->second->get_CoefficientMatchingValue("LO")*Gmu1_Gmu0;
 		
	complex_t C10_22=(0.27924*this->find("C1")->second->get_CoefficientMatchingValue("NLO")+0.33157*this->find("C4")->second->get_CoefficientMatchingValue("NLO")+2.35917*this->find("C9")->second->get_CoefficientMatchingValue("LO")
    +3.29679*this->find("C10")->second->get_CoefficientMatchingValue("LO"))*log(wilson_p("WPARAM_RUN_SM", 2))+(1.-wilson_p("WPARAM_RUN_SM", 2))*(0.26087*this->find("C9")->second->get_CoefficientMatchingValue("NLO")+1.15942*this->find("C10")->second->get_CoefficientMatchingValue("NLO"))+C1022;
	for(int ie=0;ie<=7;ie++) C10_22+=pow(wilson_p("WPARAM_RUN_SM", 2),a[ie]+2.)*((e_1a[ie]/wilson_p("WPARAM_RUN_SM", 2)+e_1b[ie])*this->find("C1")->second->get_CoefficientMatchingValue("NLO")
    +(e_4a[ie]/wilson_p("WPARAM_RUN_SM", 2)+e_4b[ie])*this->find("C4")->second->get_CoefficientMatchingValue("NLO")	
    +e_1[ie]*this->find("C1")->second->get_CoefficientMatchingValue("NNLO")
    +e_2[ie]*this->find("C2")->second->get_CoefficientMatchingValue("NNLO")
    +e_3[ie]*this->find("C3")->second->get_CoefficientMatchingValue("NNLO")
    +e_4[ie]*this->find("C4")->second->get_CoefficientMatchingValue("NNLO")
    +e_5[ie]*this->find("C5")->second->get_CoefficientMatchingValue("NNLO")
    +e_6[ie]*this->find("C6")->second->get_CoefficientMatchingValue("NNLO"));

	// C0b[10]+=alpha_ew/wilson_p("WPARAM_RUN_SM", 1)*(4.*PI/wilson_p("WPARAM_RUN_SM", 1)*C10_02+C10_12)+alpha_ew/4./PI*C10_22;

    it->second->set_WilsonCoeffRun("LO", it->second->get_CoefficientMatchingValue("LO") + 
    alpha_ew/wilson_p("WPARAM_RUN_SM", 1)*(4.*PI/wilson_p("WPARAM_RUN_SM", 1)*C10_02+C10_12)+alpha_ew/4./PI*C10_22);

    this->base["LO"] = 1;

}

void BCoefficientGroup::set_base_2_LO() {
    std::vector<std::string> coeffs {"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"};
    std::vector<complex_t> coeffs_t {};
    std::vector<complex_t> coeffs_b(10);
    coeffs_t.push_back(this->at("C1")->get_CoefficientMatchingValue("LO")/2.);

    coeffs_t.push_back(-this->at("C1")->get_CoefficientMatchingValue("LO")/6.+this->at("C2")->get_CoefficientMatchingValue("LO"));
	coeffs_t.push_back((54.*this->at("C3")->get_CoefficientMatchingValue("LO")-9.*this->at("C4")->get_CoefficientMatchingValue("LO")+864.*this->at("C5")->get_CoefficientMatchingValue("LO")-144.*this->at("C6")->get_CoefficientMatchingValue("LO"))/54.); 
	coeffs_t.push_back(this->at("C4")->get_CoefficientMatchingValue("LO")/2.+8.*this->at("C6")->get_CoefficientMatchingValue("LO"));
	coeffs_t.push_back((54.*this->at("C3")->get_CoefficientMatchingValue("LO")-9.*this->at("C4")->get_CoefficientMatchingValue("LO")+216.*this->at("C5")->get_CoefficientMatchingValue("LO")-36.*this->at("C6")->get_CoefficientMatchingValue("LO"))/54.);
 	coeffs_t.push_back(this->at("C4")->get_CoefficientMatchingValue("LO")/2.+2.*this->at("C6")->get_CoefficientMatchingValue("LO"));
	coeffs_t.push_back(this->at("C7")->get_CoefficientMatchingValue("LO"));
	coeffs_t.push_back(this->at("C8")->get_CoefficientMatchingValue("LO"));
	coeffs_t.push_back(this->at("C9")->get_CoefficientMatchingValue("LO"));
	coeffs_t.push_back(this->at("C10")->get_CoefficientMatchingValue("LO"));

    complex_t C0t7= coeffs_t[6]-1./3.*coeffs_t[4]-coeffs_t[5]; 

	complex_t C0t8= coeffs_t[7]+coeffs_t[4]; 
    for (int i=0; i<8; i++) {
        for (int j=0; j<8;j++) {
            if (j<6)
		    {
			    coeffs_b[i] += rh.get_matrix().V0[i][j]*coeffs_t[j];
		    }
		    if (j==6)
		    {
		    	coeffs_b[i] += rh.get_matrix().V0[i][j]*C0t7;
		    }
		    if (j==7)
		    {
			    coeffs_b[i] += rh.get_matrix().V0[i][j]*C0t8;
		    }
        }
    }
    for (int j=0; j<8; j++) {
        coeffs_b[8] += 4.*PI/wilson_p("WPARAM_RUN_SM", 1)*(rh.get_matrix().V0[9-1][j]*coeffs_t[j]);
    }

    coeffs_b[9] = coeffs_t[9];

    for (int i=0; i<coeffs.size(); i++) {
        this->at(coeffs[i])->set_WilsonCoeffRun("LO", coeffs_b[i]);
    }
    this->base["LO"] = 2;
}

void BCoefficientGroup::set_base_1_NLO() {
    std::vector<std::string> coeffs {"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"};
    complex_t C7_eff= this->at("C7")->get_CoefficientMatchingValue("NLO")-1./3.*this->at("C3")->get_CoefficientMatchingValue("NLO")-4./9.*this->at("C4")->get_CoefficientMatchingValue("NLO")
    -20./3.*this->at("C5")->get_CoefficientMatchingValue("NLO")-80./9.*this->at("C6")->get_CoefficientMatchingValue("NLO"); 
	complex_t C8_eff= this->at("C8")->get_CoefficientMatchingValue("NLO")+this->at("C3")->get_CoefficientMatchingValue("NLO")-1./6.*this->at("C4")->get_CoefficientMatchingValue("NLO")
    +20.*this->at("C5")->get_CoefficientMatchingValue("NLO")-10./3.*this->at("C6")->get_CoefficientMatchingValue("NLO"); 

	complex_t C7_eff_0= this->at("C7")->get_CoefficientMatchingValue("LO")-1./3.*this->at("C3")->get_CoefficientMatchingValue("LO")-4./9.*this->at("C4")->get_CoefficientMatchingValue("LO")
    -20./3.*this->at("C5")->get_CoefficientMatchingValue("LO")-80./9.*this->at("C6")->get_CoefficientMatchingValue("LO"); 
	complex_t C8_eff_0= this->at("C8")->get_CoefficientMatchingValue("LO")+this->at("C3")->get_CoefficientMatchingValue("LO")-1./6.*this->at("C4")->get_CoefficientMatchingValue("LO")
    +20.*this->at("C5")->get_CoefficientMatchingValue("LO")-10./3.*this->at("C6")->get_CoefficientMatchingValue("LO"); 

    auto calculateC1b = [&](int ie, int je, BCoefficientGroup::iterator& iterator) {
        complex_t u0_term = (rh.get_matrix().U0)[ie][je] * (je < 6 ? this->find(coeffs[je])->second->get_CoefficientMatchingValue("NLO") : (je == 6 ? C7_eff : C8_eff));
        complex_t u1_term = (rh.get_matrix().U1)[ie][je] * (je < 6 ? this->find(coeffs[je])->second->get_CoefficientMatchingValue("LO") : (je == 6 ? C7_eff_0 : C8_eff_0));
        return wilson_p("WPARAM_RUN_SM", 2) * (u0_term + u1_term);
    };


    BCoefficientGroup::iterator it = this->begin();
	for (int ie = 0; ie < 8; ie++) {
        BCoefficientGroup::iterator ite = this->begin();
        complex_t _{};
        for (int je = 0; je < 8; je++) {
            _+= calculateC1b(ie, je,  ite);
        }
        this->find(coeffs[ie])->second->set_WilsonCoeffRun("NLO", _);
    }

	double fourPiOverAlphasMu = 4.0 * PI / wilson_p("WPARAM_RUN_SM", 1);

    auto updateC1b = [&](int je, BCoefficientGroup::iterator& iterator) {
        return wilson_p("WPARAM_RUN_SM", 2) * ((rh.get_matrix().U0)[8][je] * this->find(coeffs[je])->second->get_CoefficientMatchingValue("NLO") + (rh.get_matrix().U1)[8][je] * this->find(coeffs[je])->second->get_CoefficientMatchingValue("LO"));
    };

    BCoefficientGroup::iterator iterator = this->begin();
    complex_t _{};
    for (int je = 0; je < 8; je++) {

        _ += fourPiOverAlphasMu * updateC1b(je, iterator);
        iterator++;
    }

    _ += fourPiOverAlphasMu * wilson_p("WPARAM_RUN_SM", 2) * (rh.get_matrix().U0)[8][8] * this->find("C9")->second->get_CoefficientMatchingValue("LO");
    this->find("C9")->second->set_WilsonCoeffRun("NLO", _);
    this->find("C10")->second->set_WilsonCoeffRun("NLO", wilson_p("WPARAM_RUN_SM", 2) * this->find("C10")->second->get_CoefficientMatchingValue("NLO"));
    this->base["NLO"] = 1;
}

void BCoefficientGroup::set_base_2_NLO() {

    std::vector<std::string> coeffs {"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"};
    std::vector<complex_t> coeffs_t {};
    std::vector<complex_t> coeffs_t_0 {};
    // std::vector<complex_t> coeffs_b_0 {};
    std::vector<complex_t> coeffs_b(10);

    // std::cout << "WOOOOWOOOO" << std::endl;
    coeffs_t_0.push_back(this->at("C1")->get_CoefficientMatchingValue("LO")/2.);
	coeffs_t.push_back((-5.*this->at("C1")->get_CoefficientMatchingValue("LO")+3.*(-4.*this->at("C2")->get_CoefficientMatchingValue("LO")+this->at("C1")->get_CoefficientMatchingValue("NLO")))/6.);
	coeffs_t_0.push_back(-this->at("C1")->get_CoefficientMatchingValue("LO")/6.+this->at("C2")->get_CoefficientMatchingValue("LO"));
	coeffs_t.push_back(-((11.*this->at("C1")->get_CoefficientMatchingValue("LO"))/18.)+(2.*this->at("C2")->get_CoefficientMatchingValue("LO"))/3.
    -this->at("C1")->get_CoefficientMatchingValue("NLO")/6.+this->at("C2")->get_CoefficientMatchingValue("NLO"));
	coeffs_t_0.push_back((54.*this->at("C3")->get_CoefficientMatchingValue("LO")-9.*this->at("C4")->get_CoefficientMatchingValue("LO")+864.*this->at("C5")->get_CoefficientMatchingValue("LO")
    -144.*this->at("C6")->get_CoefficientMatchingValue("LO"))/54.); 
	coeffs_t.push_back((36.*this->at("C3")->get_CoefficientMatchingValue("LO")-33.*this->at("C4")->get_CoefficientMatchingValue("LO")+4416.*this->at("C5")->get_CoefficientMatchingValue("LO")
    -700.*this->at("C6")->get_CoefficientMatchingValue("LO")+54.*this->at("C3")->get_CoefficientMatchingValue("NLO")-9.*this->at("C4")->get_CoefficientMatchingValue("NLO")
    +864.*this->at("C5")->get_CoefficientMatchingValue("NLO")-144.*this->at("C6")->get_CoefficientMatchingValue("NLO"))/54.);
	coeffs_t_0.push_back(this->at("C4")->get_CoefficientMatchingValue("LO")/2.+8.*this->at("C6")->get_CoefficientMatchingValue("LO"));
	coeffs_t.push_back((-36.*this->at("C3")->get_CoefficientMatchingValue("LO")-15.*this->at("C4")->get_CoefficientMatchingValue("LO")-960.*this->at("C5")->get_CoefficientMatchingValue("LO")
    -644.*this->at("C6")->get_CoefficientMatchingValue("LO")+9.*this->at("C4")->get_CoefficientMatchingValue("NLO")+144.*this->at("C6")->get_CoefficientMatchingValue("NLO"))/18.);
	coeffs_t_0.push_back((54.*this->at("C3")->get_CoefficientMatchingValue("LO")-9.*this->at("C4")->get_CoefficientMatchingValue("LO")+216.*this->at("C5")->get_CoefficientMatchingValue("LO")
    -36.*this->at("C6")->get_CoefficientMatchingValue("LO"))/54.);
   	coeffs_t.push_back((-36.*this->at("C3")->get_CoefficientMatchingValue("LO")+33.*this->at("C4")->get_CoefficientMatchingValue("LO")-4080.*this->at("C5")->get_CoefficientMatchingValue("LO")
    -40.*this->at("C6")->get_CoefficientMatchingValue("LO")+54.*this->at("C3")->get_CoefficientMatchingValue("NLO")-9.*this->at("C4")->get_CoefficientMatchingValue("NLO")
    +216.*this->at("C5")->get_CoefficientMatchingValue("NLO")-36.*this->at("C6")->get_CoefficientMatchingValue("NLO"))/54.);
 	coeffs_t_0.push_back(this->at("C4")->get_CoefficientMatchingValue("LO")/2.+2.*this->at("C6")->get_CoefficientMatchingValue("LO"));
 	coeffs_t.push_back((36.*this->at("C3")->get_CoefficientMatchingValue("LO")+15.*this->at("C4")->get_CoefficientMatchingValue("LO")+624.*this->at("C5")->get_CoefficientMatchingValue("LO")
    +808.*this->at("C6")->get_CoefficientMatchingValue("LO")+9.*this->at("C4")->get_CoefficientMatchingValue("NLO")+36.*this->at("C6")->get_CoefficientMatchingValue("NLO"))/18.);
	coeffs_t_0.push_back(this->at("C7")->get_CoefficientMatchingValue("LO"));
	coeffs_t.push_back(this->at("C7")->get_CoefficientMatchingValue("NLO"));
	coeffs_t_0.push_back(this->at("C8")->get_CoefficientMatchingValue("LO"));
	coeffs_t.push_back(this->at("C8")->get_CoefficientMatchingValue("NLO"));
	coeffs_t_0.push_back(this->at("C9")->get_CoefficientMatchingValue("LO"));
	coeffs_t.push_back(this->at("C9")->get_CoefficientMatchingValue("NLO"));
	coeffs_t_0.push_back(this->at("C10")->get_CoefficientMatchingValue("LO"));
	coeffs_t.push_back(this->at("C10")->get_CoefficientMatchingValue("NLO"));

    complex_t C0t7= coeffs_t_0[6]-1./3.*coeffs_t_0[4]-coeffs_t_0[5]; 
	complex_t C0t8= coeffs_t_0[7]+coeffs_t_0[4];

	complex_t C1t7= coeffs_t[6]-1./3.*coeffs_t[4]-coeffs_t[5]; 
	complex_t C1t8= coeffs_t[7]+coeffs_t[4];

    for (int i=0; i<8; i++) {
        for (int j=0; j<8;j++) {
            if (j<6)
		    {
                coeffs_b[i] += wilson_p("WPARAM_RUN_SM", 2)*(rh.get_matrix().V0[i][j]*coeffs_t[j]+rh.get_matrix().V1[i][j]*coeffs_t_0[j]);
		    }
		    if (j==6)
		    {
                coeffs_b[i] += wilson_p("WPARAM_RUN_SM", 2)*(rh.get_matrix().V0[i][j]*C1t7+rh.get_matrix().V1[i][j]*C0t7);
		    }
		    if (j==7)
		    {
                coeffs_b[i] += wilson_p("WPARAM_RUN_SM", 2)*(rh.get_matrix().V0[i][j]*C1t8+rh.get_matrix().V1[i][j]*C0t8);
		    }
        }
    }

    for (int j=0; j<8; j++) {
        coeffs_b[8] += 4.*PI/wilson_p("WPARAM_RUN_SM", 1)*(wilson_p("WPARAM_RUN_SM", 2)*(rh.get_matrix().V0[9-1][j]*coeffs_t[j]+rh.get_matrix().V1[9-1][j]*coeffs_t_0[j]));
    }

	coeffs_b[8] += 4.*PI/wilson_p("WPARAM_RUN_SM", 1)*(wilson_p("WPARAM_RUN_SM", 2)*(rh.get_matrix().V0[9-1][9-1]*coeffs_t_0[8]));

	coeffs_b[9]=wilson_p("WPARAM_RUN_SM", 2)*coeffs_t[9];

    for (int i=0; i<coeffs.size(); i++) {
        this->at(coeffs[i])->set_WilsonCoeffRun("NLO", coeffs_b[i]);
    }
    
    this->base["NLO"] = 2;
}

void BCoefficientGroup::set_base_1_NNLO() {
    std::vector<std::string> coeffs {"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"};
    complex_t C7_eff= this->at("C7")->get_CoefficientMatchingValue("NNLO")-1./3.*this->at("C3")->get_CoefficientMatchingValue("NNLO")-4./9.*this->at("C4")->get_CoefficientMatchingValue("NNLO")
    -20./3.*this->at("C5")->get_CoefficientMatchingValue("NNLO")-80./9.*this->at("C6")->get_CoefficientMatchingValue("NNLO"); 
	complex_t C8_eff= this->at("C8")->get_CoefficientMatchingValue("NNLO")+this->at("C3")->get_CoefficientMatchingValue("NNLO")-1./6.*this->at("C4")->get_CoefficientMatchingValue("NNLO")
    +20.*this->at("C5")->get_CoefficientMatchingValue("NNLO")-10./3.*this->at("C6")->get_CoefficientMatchingValue("NNLO"); 

	complex_t C7_eff_0= this->at("C7")->get_CoefficientMatchingValue("LO")-1./3.*this->at("C3")->get_CoefficientMatchingValue("LO")-4./9.*this->at("C4")->get_CoefficientMatchingValue("LO")
    -20./3.*this->at("C5")->get_CoefficientMatchingValue("LO")-80./9.*this->at("C6")->get_CoefficientMatchingValue("LO"); 
	complex_t C8_eff_0= this->at("C8")->get_CoefficientMatchingValue("LO")+this->at("C3")->get_CoefficientMatchingValue("LO")-1./6.*this->at("C4")->get_CoefficientMatchingValue("LO")
    +20.*this->at("C5")->get_CoefficientMatchingValue("LO")-10./3.*this->at("C6")->get_CoefficientMatchingValue("LO"); 

	complex_t C7_eff_1= this->at("C7")->get_CoefficientMatchingValue("NLO")-1./3.*this->at("C3")->get_CoefficientMatchingValue("NLO")-4./9.*this->at("C4")->get_CoefficientMatchingValue("NLO")
    -20./3.*this->at("C5")->get_CoefficientMatchingValue("NLO")-80./9.*this->at("C6")->get_CoefficientMatchingValue("NLO"); 
	complex_t C8_eff_1= this->at("C8")->get_CoefficientMatchingValue("NLO")+this->at("C3")->get_CoefficientMatchingValue("NLO")-1./6.*this->at("C4")->get_CoefficientMatchingValue("NLO")
    +20.*this->at("C5")->get_CoefficientMatchingValue("NLO")-10./3.*this->at("C6")->get_CoefficientMatchingValue("NLO"); 

    auto calculateC2b = [&](int ie, int je, BCoefficientGroup::iterator& iterator) {
        complex_t u0_term = (rh.get_matrix().U0)[ie][je] * (je < 6 ? this->find(coeffs[je])->second->get_CoefficientMatchingValue("NNLO") : (je == 6 ? C7_eff : C8_eff));
        complex_t u1_term = (rh.get_matrix().U1)[ie][je] * (je < 6 ? this->find(coeffs[je])->second->get_CoefficientMatchingValue("NLO") : (je == 6 ? C7_eff_1 : C8_eff_1));
        complex_t u2_term = (rh.get_matrix().U2)[ie][je] * (je < 6 ? this->find(coeffs[je])->second->get_CoefficientMatchingValue("LO") : (je == 6 ? C7_eff_0 : C8_eff_0));
        return wilson_p("WPARAM_RUN_SM", 2) * wilson_p("WPARAM_RUN_SM", 2) * (u0_term + u1_term + u2_term);
    };

    BCoefficientGroup::iterator it = this->begin();
	for (int ie = 0; ie < 8; ie++) {
        BCoefficientGroup::iterator iterator = this->begin();
        complex_t _{};
        for (int je = 0; je < 8; je++) {
            _ += calculateC2b(ie, je, iterator);
        }
        this->find(coeffs[ie])->second->set_WilsonCoeffRun("NNLO", _);
    }

	double fourPiOverAlphasMu = 4.0 * PI / wilson_p("WPARAM_RUN_SM", 1);

    auto updateC2b = [&](int je, BCoefficientGroup::iterator& iterator) {
        return wilson_p("WPARAM_RUN_SM", 2) * wilson_p("WPARAM_RUN_SM", 2) * ((rh.get_matrix().U0)[8][je] * this->find(coeffs[je])->second->get_CoefficientMatchingValue("NNLO") 
        + (rh.get_matrix().U1)[8][je] * this->find(coeffs[je])->second->get_CoefficientMatchingValue("NLO") + (rh.get_matrix().U2)[8][je] * this->find(coeffs[je])->second->get_CoefficientMatchingValue("LO"));
    };

    BCoefficientGroup::iterator iterator = this->begin();
    complex_t _{};
    for (int je = 0; je < 8; je++) {
        _ += fourPiOverAlphasMu * updateC2b(je, iterator);
        iterator++;
    }

    _ += fourPiOverAlphasMu * wilson_p("WPARAM_RUN_SM", 2) * wilson_p("WPARAM_RUN_SM", 2) * ((rh.get_matrix().U0)[8][8] * this->find("C9")->second->get_CoefficientMatchingValue("NLO") + (rh.get_matrix().U1)[8][8] * this->find("C9")->second->get_CoefficientMatchingValue("LO"));
    this->find("C9")->second->set_WilsonCoeffRun("NNLO", _);
    this->find("C10")->second->set_WilsonCoeffRun("NNLO",wilson_p("WPARAM_RUN_SM", 2) * wilson_p("WPARAM_RUN_SM", 2) * this->find("C10")->second->get_CoefficientMatchingValue("NNLO"));
    this->base["NNLO"] = 1;
}

void BCoefficientGroup::set_base_2_NNLO() {
    
}

void BScalarCoefficientGroup::set_base_1_LO() {
    complex_t coeff_temp= this->at("CQ1")->get_CoefficientMatchingValue("LO")* pow(wilson_p("WPARAM_RUN_SM", 2),-4./wilson_p("WPARAM_SI_SM", 5));
    this->at("CQ1")->set_WilsonCoeffRun("LO", coeff_temp);
    complex_t coeff_temp2= this->at("CQ2")->get_CoefficientMatchingValue("LO")* pow(wilson_p("WPARAM_RUN_SM", 2),-4./wilson_p("WPARAM_SI_SM", 5));
    this->at("CQ2")->set_WilsonCoeffRun("LO", coeff_temp2);

}

void BScalarCoefficientGroup::set_base_1_NLO() {
    complex_t coeff_temp= this->at("CQ1")->get_CoefficientMatchingValue("NLO")* pow(wilson_p("WPARAM_RUN_SM", 2),-4./wilson_p("WPARAM_SI_SM", 5))*wilson_p("WPARAM_RUN_SM", 2);
    this->at("CQ1")->set_WilsonCoeffRun("NLO", coeff_temp);
    complex_t coeff_temp2= this->at("CQ2")->get_CoefficientMatchingValue("NLO")* pow(wilson_p("WPARAM_RUN_SM", 2),-4./wilson_p("WPARAM_SI_SM", 5))*wilson_p("WPARAM_RUN_SM", 2);
    this->at("CQ2")->set_WilsonCoeffRun("NLO", coeff_temp2);
    // std::cout << coeff_temp2 << std::endl;
}

void BPrimeCoefficientGroup::set_base_1_LO() {
    complex_t coeff_temp= this->at("CP7")->get_CoefficientMatchingValue("LO")* std::pow(wilson_p("WPARAM_RUN_SM", 2), 16. / 23.) ;
    this->at("CP7")->set_WilsonCoeffRun("LO", coeff_temp);
    complex_t coeff_temp2= this->at("CP8")->get_CoefficientMatchingValue("LO")* std::pow(wilson_p("WPARAM_RUN_SM", 2), 14. / 23.);
    this->at("CP8")->set_WilsonCoeffRun("LO", coeff_temp2);

    complex_t coeff_temp5= this->at("CP9")->get_CoefficientMatchingValue("LO");
    this->at("CP9")->set_WilsonCoeffRun("LO", coeff_temp5);
    complex_t coeff_temp6= this->at("CP10")->get_CoefficientMatchingValue("LO");
    this->at("CP10")->set_WilsonCoeffRun("LO", coeff_temp6);

    complex_t coeff_temp3= this->at("CPQ1")->get_CoefficientMatchingValue("LO")* pow(wilson_p("WPARAM_RUN_SM", 2),-4./wilson_p("WPARAM_SI_SM", 5));
    this->at("CPQ1")->set_WilsonCoeffRun("LO", coeff_temp3);
    complex_t coeff_temp4= this->at("CPQ2")->get_CoefficientMatchingValue("LO")* pow(wilson_p("WPARAM_RUN_SM", 2),-4./wilson_p("WPARAM_SI_SM", 5));
    this->at("CPQ2")->set_WilsonCoeffRun("LO", coeff_temp4);

}

void CoefficientGroup::claim_coefficients() {
    for (auto& [_, coeff]: *this) {
        coeff->set_owned(true);
    }
}

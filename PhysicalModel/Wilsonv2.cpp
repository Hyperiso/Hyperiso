#include "Wilsonv2.h"
#include <iostream>


std::complex<double> C1::NLO_calculation() {
    double L=log(this->get_Q_match()*this->get_Q_match()/(*sm)("MASS",24)/(*sm)("MASS",24));
    double coeff_temp = 15.+6.*L;
    return this->double_to_complex_save("NLO", coeff_temp);
}

std::complex<double> C1::NNLO_calculation() {
    double L=log(this->get_Q_match()*this->get_Q_match()/(*sm)("MASS",24)/(*sm)("MASS",24));
    double coeff_temp = -T(W_param->xt)+7987./72.+17.*PI*PI/3.+475./6.*L+17.*L*L;
    return this->double_to_complex_save("NNLO", coeff_temp);
}



std::complex<double> C2::LO_calculation() {
    return this->double_to_complex_save("LO", 1.);
}


std::complex<double> C2::NNLO_calculation() {
    double L=log(this->get_Q_match()*this->get_Q_match()/(*sm)("MASS",24)/(*sm)("MASS",24));
    double coeff_temp = 127./18.+4./3.*PI*PI+46./3.*L+4.*L*L;
    return this->double_to_complex_save("NNLO", coeff_temp);
}

std::complex<double> C4::NLO_calculation() {
    double L=log(this->get_Q_match()*this->get_Q_match()/(*sm)("MASS",24)/(*sm)("MASS",24));
    double coeff_temp = E0t(W_param->xt)-7./9.+2./3.*L;
    return this->double_to_complex_save("NLO", coeff_temp);
}

std::complex<double> C4::NNLO_calculation() {
    double L=log(this->get_Q_match()*this->get_Q_match()/(*sm)("MASS",24)/(*sm)("MASS",24));
    double coeff_temp = E1t(W_param->xt,log(this->get_Q_match()*this->get_Q_match()/W_param->mass_top_muW/W_param->mass_top_muW))+950./243.+10./81.*PI*PI+124./27.*L+10./27.*L*L;
    return this->double_to_complex_save("NNLO", coeff_temp);
}


std::complex<double> C7::LO_calculation() {
    double coeff_temp = -0.5*A0t(W_param->xt)-23./36.;
    return this->double_to_complex_save("LO", coeff_temp);
}

std::complex<double> C7::NLO_calculation() {
    double L=log(this->get_Q_match()*this->get_Q_match()/(*sm)("MASS",24)/(*sm)("MASS",24));
    double coeff_temp = -0.5*A1t(W_param->xt,log(this->get_Q_match()*this->get_Q_match()/W_param->mass_top_muW/W_param->mass_top_muW))+713./243.+4./81.*L-4./9.*(E0t(W_param->xt)-7./9.+2./3.*L);;
    return this->double_to_complex_save("NLO", coeff_temp);
}

std::complex<double> C7::NNLO_calculation() {

    double xtW=pow(sm->running_mass((*sm)("MASS",6),(*sm)("MASS",6), (*sm)("MASS",24))/(*sm)("MASS", 24), 2); // mass top at pole for mtot param
	double xtt=pow((*sm)("MASS",6)/(*sm)("MASS",24),2.); // 24 -> W

    double coeff_temp = (C7t2mt(xtt)+log(this->get_Q_match()*this->get_Q_match()/W_param->mass_top_muW/W_param->mass_top_muW)*((-592.*pow(W_param->xt,5.)-22.*pow(W_param->xt,4.)+12814.*pow(W_param->xt,3.)-6376.*W_param->xt*W_param->xt+512.*W_param->xt)/27./pow(W_param->xt-1.,5.)*Li2(1.-1./W_param->xt)
	+(-26838.*pow(W_param->xt,5.)+25938.*pow(W_param->xt,4.)+627367.*pow(W_param->xt,3.)-331956.*W_param->xt*W_param->xt+16989.*W_param->xt-460.)/729./pow(W_param->xt-1.,6.)*log(W_param->xt)
	+(34400.*pow(W_param->xt,5.)+276644.*pow(W_param->xt,4.)-2668324.*pow(W_param->xt,3.)+1694437.*W_param->xt*W_param->xt-323354.*W_param->xt+53077.)/2187./pow(W_param->xt-1.,5.)
	+log(this->get_Q_match()*this->get_Q_match()/W_param->mass_top_muW/W_param->mass_top_muW)*((-63.*pow(W_param->xt,5.)+532.*pow(W_param->xt,4.)+2089.*pow(W_param->xt,3.)-1118.*W_param->xt*W_param->xt)/9./pow(W_param->xt-1.,6.)*log(W_param->xt)
	+(1186.*pow(W_param->xt,5.)-2705.*pow(W_param->xt,4.)-24791.*pow(W_param->xt,3.)-16099.*W_param->xt*W_param->xt+19229.*W_param->xt-2740.)/162./pow(W_param->xt-1.,5.))) )
	-(C7c2MW(xtW)+13763./2187.*log(this->get_Q_match()*this->get_Q_match()/(*sm)("MASS",24)/(*sm)("MASS",24))+814./729.*pow(log(this->get_Q_match()*this->get_Q_match()/(*sm)("MASS",24)/(*sm)("MASS",24)),2.));
    return this->double_to_complex_save("NNLO", coeff_temp);
}

std::complex<double> C8::LO_calculation() {

    double coeff_temp = -0.5*F0t(W_param->xt)-1./3.;
    return this->double_to_complex_save("LO", coeff_temp);
}

std::complex<double> C8::NLO_calculation() {
    double L=log(this->get_Q_match()*this->get_Q_match()/(*sm)("MASS",24)/(*sm)("MASS",24));
    double coeff_temp = -0.5*F1t(W_param->xt,log(this->get_Q_match()*this->get_Q_match()/W_param->mass_top_muW/W_param->mass_top_muW))+91./324.-4./27.*L-(E0t(W_param->xt)-7./9.+2./3.*L)/6.;
    return this->double_to_complex_save("NLO", coeff_temp);
}

std::complex<double> C8::NNLO_calculation() {

    double xtW=pow(sm->running_mass((*sm)("MASS",6),(*sm)("MASS",6), (*sm)("MASS",24))/(*sm)("MASS", 24), 2); // mass top at pole for mtot param
	double xtt=pow((*sm)("MASS",6)/(*sm)("MASS",24),2.); // 24 -> W

    double coeff_temp = (C8t2mt(xtt)+log(this->get_Q_match()*this->get_Q_match()/W_param->mass_top_muW/W_param->mass_top_muW)*((-148.*pow(W_param->xt,5.)+1052.*pow(W_param->xt,4.)-4811.*pow(W_param->xt,3.)-3520.*W_param->xt*W_param->xt-61.*W_param->xt)/18./pow(W_param->xt-1.,5.)*Li2(1.-1./W_param->xt)
	+(-15984.*pow(W_param->xt,5.)+152379.*pow(W_param->xt,4.)-1358060.*pow(W_param->xt,3.)-1201653.*W_param->xt*W_param->xt-74190.*W_param->xt+9188.)/1944./pow(W_param->xt-1.,6.)*log(W_param->xt)
	+(109669.*pow(W_param->xt,5.)-1112675.*pow(W_param->xt,4.)+6239377.*pow(W_param->xt,3.)+8967623.*W_param->xt*W_param->xt+768722.*W_param->xt-42796.)/11664./pow(W_param->xt-1.,5.)
	+log(this->get_Q_match()*this->get_Q_match()/W_param->mass_top_muW/W_param->mass_top_muW)*((-139.*pow(W_param->xt,4.)-2938.*pow(W_param->xt,3.)-2683.*W_param->xt*W_param->xt)/12./pow(W_param->xt-1.,6.)*log(W_param->xt)
	+(1295.*pow(W_param->xt,5.)-7009.*pow(W_param->xt,4.)+29495.*pow(W_param->xt,3.)+64513.*W_param->xt*W_param->xt+17458.*W_param->xt-2072.)/216./pow(W_param->xt-1.,5.))) )
	-(C8c2MW(xtW)+16607./5832.*log(this->get_Q_match()*this->get_Q_match()/(*sm)("MASS",24)/(*sm)("MASS",24))+397./486.*pow(log(this->get_Q_match()*this->get_Q_match()/(*sm)("MASS",24)/(*sm)("MASS",24)),2.));
    return this->double_to_complex_save("NNLO", coeff_temp);

}


std::complex<double> C9::LO_calculation() {

    double L=log(this->get_Q_match()*this->get_Q_match()/(*sm)("MASS",24)/(*sm)("MASS",24));
    double coeff_temp = (1.-4.*W_param->sw2)/W_param->sw2*C0t(W_param->xt)-B0t(W_param->xt)/W_param->sw2-D0t(W_param->xt) +1./4./W_param->sw2+38./27.-4./9.*L;
    return this->double_to_complex_save("LO", coeff_temp);
}

std::complex<double> C9::NLO_calculation() {
    double L=log(this->get_Q_match()*this->get_Q_match()/(*sm)("MASS",24)/(*sm)("MASS",24));
    double coeff_temp = (1.-4.*W_param->sw2)/W_param->sw2*C1t(W_param->xt,log(this->get_Q_match()*this->get_Q_match()/W_param->mass_top_muW/W_param->mass_top_muW))
    -B1t(W_param->xt,log(this->get_Q_match()*this->get_Q_match()/W_param->mass_top_muW/W_param->mass_top_muW))/W_param->sw2
    -D1t(W_param->xt,log(this->get_Q_match()*this->get_Q_match()/W_param->mass_top_muW/W_param->mass_top_muW)) +1./W_param->sw2+524./729.-128./243.*PI*PI-16./3.*L-128./81.*L*L;
    return this->double_to_complex_save("NLO", coeff_temp);
}

std::complex<double> C9::NNLO_calculation() {

    this->set_CoefficientMatchingValue("LO", 1);
    return 1;
}


std::complex<double> C10::LO_calculation() {

    double coeff_temp = (B0t(W_param->xt)-C0t(W_param->xt))/W_param->sw2-1./4./W_param->sw2;
    return this->double_to_complex_save("LO", coeff_temp);
}

std::complex<double> C10::NLO_calculation() {
    double coeff_temp =  (B1t(W_param->xt,log(this->get_Q_match()*this->get_Q_match()/W_param->mass_top_muW/W_param->mass_top_muW))
    -C1t(W_param->xt,log(this->get_Q_match()*this->get_Q_match()/W_param->mass_top_muW/W_param->mass_top_muW)))/W_param->sw2-1./W_param->sw2;
    return this->double_to_complex_save("NLO", coeff_temp);
}

std::complex<double> C10::NNLO_calculation() {

    this->set_CoefficientMatchingValue("LO", 1);
    return 1;
}

void BCoefficientGroup::set_base_1_LO() {
    if (this->at("C1")->get_Q() == 0.) {
        LOG_ERROR("BORDELDESALOPERIEDECODEDEMERDE", "IL FAUT METTRE Q AVANT DE LANCER CETTE *****");
    }
    complex_t C7_eff= this->at("C7")->get_CoefficientMatchingValue("LO")-1./3.*this->at("C3")->get_CoefficientMatchingValue("LO")
    -4./9.*this->at("C4")->get_CoefficientMatchingValue("LO")-20./3.*this->at("C5")->get_CoefficientMatchingValue("LO")-80./9.*this->at("C6")->get_CoefficientMatchingValue("LO"); 
	complex_t C8_eff= this->at("C8")->get_CoefficientMatchingValue("LO")+this->at("C3")->get_CoefficientMatchingValue("LO")-1./6.*this->at("C4")->get_CoefficientMatchingValue("LO")
    +20.*this->at("C5")->get_CoefficientMatchingValue("LO")-10./3.*this->at("C6")->get_CoefficientMatchingValue("LO"); 

    std::cout << C7_eff << std::endl;
    std::cout << C8_eff << std::endl;

    auto calculateC0b = [&](int ie, int je, BCoefficientGroup::iterator& iterator) {
        return (W_param->U0)[ie][je] * (je < 6 ? (iterator++)->second->get_CoefficientMatchingValue("LO") : (je == 6 ? C7_eff : C8_eff));
    };

    BCoefficientGroup::iterator it = this->begin();
    for (int ie = 0; ie < 8; ie++) {
        complex_t _{};
        BCoefficientGroup::iterator it2 = this->begin();
        for (int je = 0; je < 8; je++) {
            _+= calculateC0b(ie, je, it2);

        }
        it->second->set_WilsonCoeffRun("LO", _);
        it++;
    }
    
	double fourPiOverAlphasMu = 4.0 * PI / W_param->alphas_mu;

    auto updateC0b = [&](int je, BCoefficientGroup::iterator& iterator) {
        return (W_param->U0)[8][je] * (iterator++)->second->get_CoefficientMatchingValue("LO");
    };
    
    {
        complex_t _{};
        BCoefficientGroup::iterator itera = this->begin(); 
        for (int je = 0; je < 8; je++) {
        _ += fourPiOverAlphasMu * updateC0b(je, itera);
        }
        (it++)->second->set_WilsonCoeffRun("LO", _);
    }

    it->second->set_WilsonCoeffRun("LO", it->second->get_CoefficientMatchingValue("LO"));

    if (++it == this->end()) {
        std::cout << "Ca marche ! " << std::endl;
    }
}

void BCoefficientGroup::set_base_1_NLO() {
    complex_t C7_eff= this->at("C7")->get_CoefficientMatchingValue("NLO")-1./3.*this->at("C3")->get_CoefficientMatchingValue("NLO")-4./9.*this->at("C4")->get_CoefficientMatchingValue("NLO")
    -20./3.*this->at("C5")->get_CoefficientMatchingValue("NLO")-80./9.*this->at("C6")->get_CoefficientMatchingValue("NLO"); 
	complex_t C8_eff= this->at("C8")->get_CoefficientMatchingValue("NLO")+this->at("C3")->get_CoefficientMatchingValue("NLO")-1./6.*this->at("C4")->get_CoefficientMatchingValue("NLO")
    +20.*this->at("C5")->get_CoefficientMatchingValue("NLO")-10./3.*this->at("C6")->get_CoefficientMatchingValue("NLO"); 

	complex_t C7_eff_0= this->at("C7")->get_CoefficientMatchingValue("LO")-1./3.*this->at("C3")->get_CoefficientMatchingValue("LO")-4./9.*this->at("C4")->get_CoefficientMatchingValue("LO")
    -20./3.*this->at("C5")->get_CoefficientMatchingValue("LO")-80./9.*this->at("C6")->get_CoefficientMatchingValue("LO"); 
	complex_t C8_eff_0= this->at("C8")->get_CoefficientMatchingValue("LO")+this->at("C3")->get_CoefficientMatchingValue("LO")-1./6.*this->at("C4")->get_CoefficientMatchingValue("LO")
    +20.*this->at("C5")->get_CoefficientMatchingValue("LO")-10./3.*this->at("C6")->get_CoefficientMatchingValue("LO"); 

    auto calculateC1b = [&](int ie, int je, BCoefficientGroup::iterator& iterator) {
        complex_t u0_term = (W_param->U0)[ie][je] * (je < 6 ? iterator->second->get_CoefficientMatchingValue("NLO") : (je == 6 ? C7_eff : C8_eff));
        complex_t u1_term = (W_param->U1)[ie][je] * (je < 6 ? (iterator++)->second->get_CoefficientMatchingValue("LO") : (je == 6 ? C7_eff_0 : C8_eff_0));
        return W_param->eta_mu * (u0_term + u1_term);
    };


    BCoefficientGroup::iterator it = this->begin();
	for (int ie = 0; ie < 8; ie++) {
        BCoefficientGroup::iterator ite = this->begin();
        complex_t _{};
        for (int je = 0; je < 8; je++) {
            _+= calculateC1b(ie, je,  ite);
        }
        (it++)->second->set_WilsonCoeffRun("NLO", _);
    }

	double fourPiOverAlphasMu = 4.0 * PI / W_param->alphas_mu;

    auto updateC1b = [&](int je, BCoefficientGroup::iterator& iterator) {
        return W_param->eta_mu * ((W_param->U0)[8][je] * iterator->second->get_CoefficientMatchingValue("NLO") + (W_param->U1)[8][je] * iterator->second->get_CoefficientMatchingValue("LO"));
    };

    BCoefficientGroup::iterator iterator = this->begin();
    complex_t _{};
    for (int je = 0; je < 8; je++) {

        _ += fourPiOverAlphasMu * updateC1b(je, iterator);
        iterator++;
    }

    _ += fourPiOverAlphasMu * W_param->eta_mu * (W_param->U0)[8][8] * it->second->get_CoefficientMatchingValue("NLO");
    (it++)->second->set_WilsonCoeffRun("NLO", _);
    it->second->set_WilsonCoeffRun("NLO", W_param->eta_mu * it->second->get_CoefficientMatchingValue("NLO"));
}

int main() {
    MemoryManager::GetInstance()->init();

    C4 C4_test{81.};
    C4_test.LO_calculation();
    C4_test.NLO_calculation();
    C4_test.NNLO_calculation();

    BCoefficientGroup bcoeff{81.};
    bcoeff.set_Q_run(81.);
    bcoeff.init_LO();
    bcoeff.set_base_1_LO();

    std::cout <<  bcoeff << std::endl;

    std::cout << ".................................................................................................." << std::endl;

    bcoeff.init_NLO();
    bcoeff.set_base_1_NLO();

    return 0;
}
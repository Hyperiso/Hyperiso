#include "Wilsonv2.h"



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

std::complex<double> C3::NNLO_calculation() {
    double L=log(this->get_Q_match()*this->get_Q_match()/(*sm)("MASS",24)/(*sm)("MASS",24));
    double coeff_temp = G1t(W_param->xt,log(this->get_Q_match()*this->get_Q_match()/W_param->mass_top_muW/W_param->mass_top_muW))-680./243.-20./81.*PI*PI-68./81.*L-20./27.*L*L;
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

std::complex<double> C5::NNLO_calculation() {
    double L=log(this->get_Q_match()*this->get_Q_match()/(*sm)("MASS",24)/(*sm)("MASS",24));
    double coeff_temp = -G1t(W_param->xt,log(this->get_Q_match()*this->get_Q_match()/W_param->mass_top_muW/W_param->mass_top_muW))/10.+2./15.*E0t(W_param->xt)+68./243.+2./81.*PI*PI+14./81.*L+2./27.*L*L;
    return this->double_to_complex_save("NNLO", coeff_temp);
}

std::complex<double> C6::NNLO_calculation() {
    double L=log(this->get_Q_match()*this->get_Q_match()/(*sm)("MASS",24)/(*sm)("MASS",24));
    double coeff_temp = -3./16.*G1t(W_param->xt,log(this->get_Q_match()*this->get_Q_match()/W_param->mass_top_muW/W_param->mass_top_muW))+E0t(W_param->xt)/4.+85./162.+5./108.*PI*PI+35./108.*L+5./36.*L*L;
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

    double xtW=pow(QCDHelper::msbar_mass(6, (*sm)("MASS",24))/(*sm)("MASS", 24), 2); // mass top at pole for mtot param
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

    double xtW=pow(QCDHelper::msbar_mass(6, (*sm)("MASS",24))/(*sm)("MASS", 24), 2); // mass top at pole for mtot param
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
    // std::cout << coeff_temp << std::endl;
    return this->double_to_complex_save("LO", coeff_temp);
}

std::complex<double> C10::NLO_calculation() {
    double coeff_temp =  (B1t(W_param->xt,log(this->get_Q_match()*this->get_Q_match()/W_param->mass_top_muW/W_param->mass_top_muW))
    -C1t(W_param->xt,log(this->get_Q_match()*this->get_Q_match()/W_param->mass_top_muW/W_param->mass_top_muW)))/W_param->sw2-1./W_param->sw2;
    return this->double_to_complex_save("NLO", coeff_temp);
}
 
std::complex<double> C10::NNLO_calculation() {
    double xtW=pow(QCDHelper::msbar_mass(6, (*sm)("MASS",24)) / (*sm)("MASS", 24), 2); // mass top at pole for mtot param
	double xtt=pow((*sm)("MASS",6)/(*sm)("MASS",24),2.); // 24 -> W
    double coeff_temp = ((C10Wt2mt(xtt)+log(this->get_Q_match()*this->get_Q_match()/W_param->mass_top_muW/W_param->mass_top_muW)*((69.+1292.*W_param->xt-209.*W_param->xt*W_param->xt)/18./pow(W_param->xt-1.,3.)
	-(521.*W_param->xt+105.*W_param->xt*W_param->xt-50.*pow(W_param->xt,3.))/9./pow(W_param->xt-1.,4.)*log(W_param->xt)
	-(47.*W_param->xt+W_param->xt*W_param->xt)/3./pow(W_param->xt-1.,3.)*Li2(1.-1./W_param->xt)
	+log(this->get_Q_match()*this->get_Q_match()/W_param->mass_top_muW/W_param->mass_top_muW)*((61.*W_param->xt+11.*W_param->xt*W_param->xt)/3./pow(W_param->xt-1.,3.)-(49.*W_param->xt+96.*W_param->xt*W_param->xt-pow(W_param->xt,3.))/6./pow(W_param->xt-1.,4.)*log(W_param->xt))))
	-(C10Wc2MW(xtW)-23./6.*log(this->get_Q_match()*this->get_Q_match()/(*sm)("MASS", 24)/(*sm)("MASS",24)))
	+(C10Zt2mt(xtt)+log(this->get_Q_match()*this->get_Q_match()/W_param->mass_top_muW/W_param->mass_top_muW)*((188.*W_param->xt+4.*W_param->xt*W_param->xt+95.*pow(W_param->xt,3.)-47.*pow(W_param->xt,4.))/6./pow(W_param->xt-1.,3.)*Li2(1.-1./W_param->xt)
	+(1468.*W_param->xt+1578.*W_param->xt*W_param->xt-25.*pow(W_param->xt,3.)-141.*pow(W_param->xt,4.))/18./pow(W_param->xt-1.,4.)*log(W_param->xt)
	-(4622.*W_param->xt+1031.*W_param->xt*W_param->xt+582.*pow(W_param->xt,3.)-475.*pow(W_param->xt,4.))/36./pow(W_param->xt-1.,3.)
	+log(this->get_Q_match()*this->get_Q_match()/W_param->mass_top_muW/W_param->mass_top_muW)*((49.*W_param->xt+315.*W_param->xt*W_param->xt-4.*pow(W_param->xt,3.))/6./pow(W_param->xt-1.,4.)*log(W_param->xt)-(440.*W_param->xt+257.*W_param->xt*W_param->xt+72.*pow(W_param->xt,3.)-49.*pow(W_param->xt,4.))/12./pow(W_param->xt-1.,3.))))
	+C10Z2tri(xtt)
	)*(-2./W_param->sw2);
    return this->double_to_complex_save("NNLO", coeff_temp);
}
std::complex<double> CQ1::LO_calculation() {
    double CSc_SM=-W_param->xt*(W_param->xt-2.)/12./(W_param->xt-1.)/(W_param->xt-1.)+(W_param->xt-2.)*(3.*W_param->xt-1.)/24./pow(W_param->xt-1.,3.)*log(W_param->xt);
    double CSn_SMonly=-3.*W_param->xt/8./W_param->xh+W_param->xt*F0SP(W_param->xt);
    double coeff_temp=(CSc_SM+CSn_SMonly)*(W_param->ml*W_param->mass_b_muW/(*sm)("MASS",24)/(*sm)("MASS",24))/W_param->sw2;
    return this->double_to_complex_save("LO", coeff_temp);
}

std::complex<double> CQ2::LO_calculation() {
    double CPc_SM=1./24.*(W_param->xt*(36.*W_param->xt3-203.*W_param->xt2+352.*W_param->xt-209.)/6./pow(W_param->xt-1.,3.)+(17.*W_param->xt4-34.*W_param->xt3+4.*W_param->xt2+23.*W_param->xt-6.)/pow(W_param->xt-1.,4.)*log(W_param->xt))
	-W_param->sw2/36.*(W_param->xt*(18.*W_param->xt3-139.*W_param->xt2+274.*W_param->xt-129.)/2./pow(W_param->xt-1.,3.)+(24.*W_param->xt4-33.*W_param->xt3-45.*W_param->xt2+50.*W_param->xt-8.)/pow(W_param->xt-1.,4.)*log(W_param->xt));
    double CPn_SMonly=0.;
    double coeff_temp = (CPc_SM+CPn_SMonly)*(W_param->ml*W_param->mass_b_muW/(*sm)("MASS",24)/(*sm)("MASS",24))/W_param->sw2;
    return this->double_to_complex_save("LO", coeff_temp);
}

std::complex<double> CP7::LO_calculation() {
    double coeff_temp = (*sm)("MASS", 3) / W_param->mass_b_muW * (-0.5 * A0t(W_param->xt) - 23. / 36.);
    return this->double_to_complex_save("LO", coeff_temp);
}

std::complex<double> CP8::LO_calculation() {
    double coeff_temp = (*sm)("MASS", 3) / W_param->mass_b_muW * (-0.5 * F0t(W_param->xt) - 1. / 3.);
    return this->double_to_complex_save("LO", coeff_temp);
}

void BCoefficientGroup::set_base_1_LO() {
    if (this->at("C1")->get_Q() == 0.) {
        LOG_ERROR("BORDELDESALOPERIEDECODEDEMERDE", "IL FAUT METTRE Q AVANT DE LANCER CETTE *****");
    }
    complex_t C7_eff= this->at("C7")->get_CoefficientMatchingValue("LO")-1./3.*this->at("C3")->get_CoefficientMatchingValue("LO")
    -4./9.*this->at("C4")->get_CoefficientMatchingValue("LO")-20./3.*this->at("C5")->get_CoefficientMatchingValue("LO")-80./9.*this->at("C6")->get_CoefficientMatchingValue("LO"); 
	complex_t C8_eff= this->at("C8")->get_CoefficientMatchingValue("LO")+this->at("C3")->get_CoefficientMatchingValue("LO")-1./6.*this->at("C4")->get_CoefficientMatchingValue("LO")
    +20.*this->at("C5")->get_CoefficientMatchingValue("LO")-10./3.*this->at("C6")->get_CoefficientMatchingValue("LO"); 

    // std::cout << C7_eff << std::endl;
    // std::cout << C8_eff << std::endl;

    auto calculateC0b = [&](int ie, int je, std::vector<std::string>& coeff_loop) {
        return (W_param->U0)[ie][je] * (je < 6 ? this->find(coeff_loop[je])->second->get_CoefficientMatchingValue("LO") : (je == 6 ? C7_eff : C8_eff));
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

	double fourPiOverAlphasMu = 4.0 * PI / W_param->alphas_mu;

    auto updateC0b = [&](int je, BCoefficientGroup::iterator& iterator) {
        // std::cout << "truc : " << je << " " << this->find(coeff_loop[je])->first << std::endl;
        return (W_param->U0)[8][je] * this->find(coeff_loop[je])->second->get_CoefficientMatchingValue("LO");
    };
    
    {
        complex_t _{};
        BCoefficientGroup::iterator itera = this->begin(); 
        for (int je = 0; je < 8; je++) {
        _ += fourPiOverAlphasMu * updateC0b(je, itera);
        }
        (++it)->second->set_WilsonCoeffRun("LO", _);
    }
    // std::cout << "after ++it first : " << it->first << std::endl;
    // std::cout << "Niels à tord : " << it->first << " " << it->second->get_CoefficientMatchingValue("LO") << std::endl;
    // std::cout << "end of truc " << this->end()->first << std::endl;
    it = this->find("C10");
    // std::cout << "after second ++it : " << it->first << std::endl;
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
	for(int ie=0;ie<=7;ie++) C10_02+=b[ie]*pow(W_param->eta_mu,a[ie])*this->find("C2")->second->get_CoefficientMatchingValue("LO");
	
	complex_t C10_12=-0.11060*log(W_param->eta_mu)/W_param->eta_mu*this->find("C2")->second->get_CoefficientMatchingValue("LO")
    +(1./W_param->eta_mu-1.)*(0.26087*this->find("C9")->second->get_CoefficientMatchingValue("LO")+1.15942*this->find("C10")->second->get_CoefficientMatchingValue("LO"));
	for(int ie=0;ie<=7;ie++) C10_12+=pow(W_param->eta_mu,a[ie]+1.)*((d_2a[ie]/W_param->eta_mu+d_2b[ie])*this->find("C2")->second->get_CoefficientMatchingValue("LO")
    +d_1[ie]*this->find("C1")->second->get_CoefficientMatchingValue("NLO")+d_4[ie]*this->find("C4")->second->get_CoefficientMatchingValue("NLO"));
		
	double Delta_alpha=0.06; 
	double Delta_rhosw2=-0.03; 
	double Delta_rem=0.01;
	double Deltar=Delta_alpha+Delta_rhosw2+Delta_rem; 
	
	double Gmu1_Gmu0=4.*PI/alpha_ew*Deltar;
	
	complex_t C1022=(46.9287715663914-3.102350691200236*log(this->get_Q_match()*this->get_Q_match())+0.0992974073578769*log(this->get_Q_match()*this->get_Q_match())*log(this->get_Q_match()*this->get_Q_match())+0.175877*(W_param->mass_top_muW-163.5)+0.0173725*(MH-125.9))/sw2OS;
 	C1022+=-this->find("C10")->second->get_CoefficientMatchingValue("LO")*Gmu1_Gmu0;
 		
	complex_t C10_22=(0.27924*this->find("C1")->second->get_CoefficientMatchingValue("NLO")+0.33157*this->find("C4")->second->get_CoefficientMatchingValue("NLO")+2.35917*this->find("C9")->second->get_CoefficientMatchingValue("LO")
    +3.29679*this->find("C10")->second->get_CoefficientMatchingValue("LO"))*log(W_param->eta_mu)+(1.-W_param->eta_mu)*(0.26087*this->find("C9")->second->get_CoefficientMatchingValue("NLO")+1.15942*this->find("C10")->second->get_CoefficientMatchingValue("NLO"))+C1022;
	for(int ie=0;ie<=7;ie++) C10_22+=pow(W_param->eta_mu,a[ie]+2.)*((e_1a[ie]/W_param->eta_mu+e_1b[ie])*this->find("C1")->second->get_CoefficientMatchingValue("NLO")
    +(e_4a[ie]/W_param->eta_mu+e_4b[ie])*this->find("C4")->second->get_CoefficientMatchingValue("NLO")	
    +e_1[ie]*this->find("C1")->second->get_CoefficientMatchingValue("NNLO")
    +e_2[ie]*this->find("C2")->second->get_CoefficientMatchingValue("NNLO")
    +e_3[ie]*this->find("C3")->second->get_CoefficientMatchingValue("NNLO")
    +e_4[ie]*this->find("C4")->second->get_CoefficientMatchingValue("NNLO")
    +e_5[ie]*this->find("C5")->second->get_CoefficientMatchingValue("NNLO")
    +e_6[ie]*this->find("C6")->second->get_CoefficientMatchingValue("NNLO"));

	// C0b[10]+=alpha_ew/W_param->alphas_mu*(4.*PI/W_param->alphas_mu*C10_02+C10_12)+alpha_ew/4./PI*C10_22;

    it->second->set_WilsonCoeffRun("LO", it->second->get_CoefficientMatchingValue("LO") + 
    alpha_ew/W_param->alphas_mu*(4.*PI/W_param->alphas_mu*C10_02+C10_12)+alpha_ew/4./PI*C10_22);

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
			    coeffs_b[i] += W_param->V0[i][j]*coeffs_t[j];
		    }
		    if (j==6)
		    {
		    	coeffs_b[i] += W_param->V0[i][j]*C0t7;
		    }
		    if (j==7)
		    {
			    coeffs_b[i] += W_param->V0[i][j]*C0t8;
		    }
        }
    }
    for (int j=0; j<8; j++) {
        coeffs_b[8] += 4.*PI/W_param->alphas_mu*(W_param->V0[9-1][j]*coeffs_t[j]);
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
        complex_t u0_term = (W_param->U0)[ie][je] * (je < 6 ? this->find(coeffs[je])->second->get_CoefficientMatchingValue("NLO") : (je == 6 ? C7_eff : C8_eff));
        complex_t u1_term = (W_param->U1)[ie][je] * (je < 6 ? this->find(coeffs[je])->second->get_CoefficientMatchingValue("LO") : (je == 6 ? C7_eff_0 : C8_eff_0));
        return W_param->eta_mu * (u0_term + u1_term);
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

	double fourPiOverAlphasMu = 4.0 * PI / W_param->alphas_mu;

    auto updateC1b = [&](int je, BCoefficientGroup::iterator& iterator) {
        return W_param->eta_mu * ((W_param->U0)[8][je] * this->find(coeffs[je])->second->get_CoefficientMatchingValue("NLO") + (W_param->U1)[8][je] * this->find(coeffs[je])->second->get_CoefficientMatchingValue("LO"));
    };

    BCoefficientGroup::iterator iterator = this->begin();
    complex_t _{};
    for (int je = 0; je < 8; je++) {

        _ += fourPiOverAlphasMu * updateC1b(je, iterator);
        iterator++;
    }

    _ += fourPiOverAlphasMu * W_param->eta_mu * (W_param->U0)[8][8] * this->find("C9")->second->get_CoefficientMatchingValue("LO");
    this->find("C9")->second->set_WilsonCoeffRun("NLO", _);
    this->find("C10")->second->set_WilsonCoeffRun("NLO", W_param->eta_mu * this->find("C10")->second->get_CoefficientMatchingValue("NLO"));
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
                coeffs_b[i] += W_param->eta_mu*(W_param->V0[i][j]*coeffs_t[j]+W_param->V1[i][j]*coeffs_t_0[j]);
		    }
		    if (j==6)
		    {
                coeffs_b[i] += W_param->eta_mu*(W_param->V0[i][j]*C1t7+W_param->V1[i][j]*C0t7);
		    }
		    if (j==7)
		    {
                coeffs_b[i] += W_param->eta_mu*(W_param->V0[i][j]*C1t8+W_param->V1[i][j]*C0t8);
		    }
        }
    }

    for (int j=0; j<8; j++) {
        coeffs_b[8] += 4.*PI/W_param->alphas_mu*(W_param->eta_mu*(W_param->V0[9-1][j]*coeffs_t[j]+W_param->V1[9-1][j]*coeffs_t_0[j]));
    }

	coeffs_b[8] += 4.*PI/W_param->alphas_mu*(W_param->eta_mu*(W_param->V0[9-1][9-1]*coeffs_t_0[8]));

	coeffs_b[9]=W_param->eta_mu*coeffs_t[9];

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
        complex_t u0_term = (W_param->U0)[ie][je] * (je < 6 ? this->find(coeffs[je])->second->get_CoefficientMatchingValue("NNLO") : (je == 6 ? C7_eff : C8_eff));
        complex_t u1_term = (W_param->U1)[ie][je] * (je < 6 ? this->find(coeffs[je])->second->get_CoefficientMatchingValue("NLO") : (je == 6 ? C7_eff_1 : C8_eff_1));
        complex_t u2_term = (W_param->U2)[ie][je] * (je < 6 ? this->find(coeffs[je])->second->get_CoefficientMatchingValue("LO") : (je == 6 ? C7_eff_0 : C8_eff_0));
        return W_param->eta_mu * W_param->eta_mu * (u0_term + u1_term + u2_term);
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

	double fourPiOverAlphasMu = 4.0 * PI / W_param->alphas_mu;

    auto updateC2b = [&](int je, BCoefficientGroup::iterator& iterator) {
        return W_param->eta_mu * W_param->eta_mu * ((W_param->U0)[8][je] * this->find(coeffs[je])->second->get_CoefficientMatchingValue("NNLO") 
        + (W_param->U1)[8][je] * this->find(coeffs[je])->second->get_CoefficientMatchingValue("NLO") + (W_param->U2)[8][je] * this->find(coeffs[je])->second->get_CoefficientMatchingValue("LO"));
    };

    BCoefficientGroup::iterator iterator = this->begin();
    complex_t _{};
    for (int je = 0; je < 8; je++) {
        _ += fourPiOverAlphasMu * updateC2b(je, iterator);
        iterator++;
    }

    _ += fourPiOverAlphasMu * W_param->eta_mu * W_param->eta_mu * ((W_param->U0)[8][8] * this->find("C9")->second->get_CoefficientMatchingValue("NLO") + (W_param->U1)[8][8] * this->find("C9")->second->get_CoefficientMatchingValue("LO"));
    this->find("C9")->second->set_WilsonCoeffRun("NNLO", _);
    this->find("C10")->second->set_WilsonCoeffRun("NNLO",W_param->eta_mu * W_param->eta_mu * this->find("C10")->second->get_CoefficientMatchingValue("NNLO"));
    this->base["NNLO"] = 1;
}

void BCoefficientGroup::set_base_2_NNLO() {
    
}

void BScalarCoefficientGroup::set_base_1_LO() {
    complex_t coeff_temp= this->at("CQ1")->get_CoefficientMatchingValue("LO")* pow(W_param->eta_mu,-4./W_param->beta0);
    this->at("CQ1")->set_WilsonCoeffRun("LO", coeff_temp);
    complex_t coeff_temp2= this->at("CQ2")->get_CoefficientMatchingValue("LO")* pow(W_param->eta_mu,-4./W_param->beta0);
    this->at("CQ2")->set_WilsonCoeffRun("LO", coeff_temp2);

}

void BScalarCoefficientGroup::set_base_1_NLO() {
    complex_t coeff_temp= this->at("CQ1")->get_CoefficientMatchingValue("NLO")* pow(W_param->eta_mu,-4./W_param->beta0)*W_param->eta_mu;
    this->at("CQ1")->set_WilsonCoeffRun("NLO", coeff_temp);
    complex_t coeff_temp2= this->at("CQ2")->get_CoefficientMatchingValue("NLO")* pow(W_param->eta_mu,-4./W_param->beta0)*W_param->eta_mu;
    this->at("CQ2")->set_WilsonCoeffRun("NLO", coeff_temp2);
    // std::cout << coeff_temp2 << std::endl;
}

void BPrimeCoefficientGroup::set_base_1_LO() {
    complex_t coeff_temp= this->at("CP7")->get_CoefficientMatchingValue("LO")* std::pow(W_param->eta_mu, 16. / 23.) ;
    this->at("CP7")->set_WilsonCoeffRun("LO", coeff_temp);
    complex_t coeff_temp2= this->at("CP8")->get_CoefficientMatchingValue("LO")* std::pow(W_param->eta_mu, 14. / 23.);
    this->at("CP8")->set_WilsonCoeffRun("LO", coeff_temp2);

    complex_t coeff_temp5= this->at("CP9")->get_CoefficientMatchingValue("LO");
    this->at("CP9")->set_WilsonCoeffRun("LO", coeff_temp5);
    complex_t coeff_temp6= this->at("CP10")->get_CoefficientMatchingValue("LO");
    this->at("CP10")->set_WilsonCoeffRun("LO", coeff_temp6);

    complex_t coeff_temp3= this->at("CPQ1")->get_CoefficientMatchingValue("LO")* pow(W_param->eta_mu,-4./W_param->beta0);
    this->at("CPQ1")->set_WilsonCoeffRun("LO", coeff_temp3);
    complex_t coeff_temp4= this->at("CPQ2")->get_CoefficientMatchingValue("LO")* pow(W_param->eta_mu,-4./W_param->beta0);
    this->at("CPQ2")->set_WilsonCoeffRun("LO", coeff_temp4);

}

bool WilsonCoefficient::fill_from_flha() {
    if (!from_lha && MemoryManager::GetInstance()->hasWilsons()) {
        auto wc = Parameters::GetInstance(ParameterType::WILSON);
        BWilsonCoefficients id = WCoefMapper::enum_elt(this->get_name());
        Model m = MemoryManager::GetInstance()->getModel();
        int w_type = (*wc)("REWCOEF", -2);
        if (m != Model::SM && w_type == 0) {
            LOG_ERROR("Value", "SM Wilsons coefficients were given, but the selected model is not SM.");
        }

        this->set_CoefficientMatchingValue("LO", complex_t((*wc)("REWCOEF", (int)id * 10 + (int)QCDOrder::LO-1), 
                                                            (*wc)("IMWCOEF", (int)id * 10 + (int)QCDOrder::LO-1)));
        this->set_CoefficientMatchingValue("NLO", complex_t((*wc)("REWCOEF", (int)id * 10 + (int)QCDOrder::NLO-1), 
                                                             (*wc)("IMWCOEF", (int)id * 10 + (int)QCDOrder::NLO-1)));
        this->set_CoefficientMatchingValue("NNLO", complex_t((*wc)("REWCOEF", (int)id * 10 + (int)QCDOrder::NNLO-1), 
                                                              (*wc)("IMWCOEF", (int)id * 10 + (int)QCDOrder::NNLO-1)));         
        from_lha = true;
    }
    return from_lha;
}

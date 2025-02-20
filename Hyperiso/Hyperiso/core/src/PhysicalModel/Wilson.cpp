#include "Wilson.h"



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
	double xtt=pow(QCDHelper::mass_t_msbar()/(*sm)("MASS",24),2.); // 24 -> W

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
	double xtt=pow(QCDHelper::mass_t_msbar()/(*sm)("MASS",24),2.); // 24 -> W

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
	double xtt=pow(QCDHelper::mass_t_msbar()/(*sm)("MASS",24),2.); // 24 -> W
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


bool WilsonCoefficient::fill_from_flha() {
    if (!from_lha && MemoryManager::GetInstance()->hasWilsons()) {
        auto wc = Parameters::GetInstance(ParameterType::WILSON);
        WCoef id = WCoefMapper::enum_elt(this->get_name());
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

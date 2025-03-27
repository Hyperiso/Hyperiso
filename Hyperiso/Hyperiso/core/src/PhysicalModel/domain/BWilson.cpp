#include "BWilson.h"

void C1::NLO_calculation() {
    std::unordered_set<ParamId> sources {
        {ParameterType::WILSON, "WPARAM_MATCH_SM", 3}  // L
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        auto L = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 3})->get_val();
        dep_param->set_expected(15 + 6 * L);
    };

    WilsonParamComposer().compose_parameter(ParamId{ParameterType::WILSON, "B_MATCH", LhaID(3040405, 6161, 0, 0)}, sources, func);
}

void C1::NNLO_calculation() {
    std::unordered_set<ParamId> sources {
        {ParameterType::WILSON, "WPARAM_MATCH_SM", 3},              // L
        {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2, 1)}     // x_t
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        auto L = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 3})->get_val();
        auto xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
        dep_param->set_expected(-T(xt) + 7987./72. + 17. * PI2 / 3. + 475./6. * L + 17. * L * L);
    };

    WilsonParamComposer().compose_parameter(ParamId{ParameterType::WILSON, "B_MATCH", LhaID(3040405, 6161, 1, 0)}, sources, func);
}

void C2::LO_calculation() {
    // return this->double_to_complex_save("LO", 1.);
    std::unordered_set<ParamId> sources {
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        dep_param->set_expected(1.);
    };

    WilsonParamComposer().compose_parameter(ParamId{ParameterType::WILSON, "B_MATCH", LhaID(3040405, 4141, 0, 0)}, sources, func);


}

void C2::NNLO_calculation() {
    double L = wilson_p("WPARAM_MATCH_SM", 3);
    double coeff_temp = 127./18.+4./3.*PI*PI+46./3.*L+4.*L*L;
}

void C3::NNLO_calculation() {
    double L = wilson_p("WPARAM_MATCH_SM", 3);
    double coeff_temp = G1t(wilson_p("WPARAM_MATCH_SM", {2,1}),log(this->get_Q_match()*this->get_Q_match()/wilson_p("WPARAM_MATCH_SM", 6)/wilson_p("WPARAM_MATCH_SM", 6)))-680./243.-20./81.*PI*PI-68./81.*L-20./27.*L*L;
}
void C4::NLO_calculation() {
    double L = wilson_p("WPARAM_MATCH_SM", 3);
    double coeff_temp = E0t(wilson_p("WPARAM_MATCH_SM", {2,1}))-7./9.+2./3.*L;
}

void C4::NNLO_calculation() {
    double L = wilson_p("WPARAM_MATCH_SM", 3);
    double coeff_temp = E1t(wilson_p("WPARAM_MATCH_SM", {2,1}),log(this->get_Q_match()*this->get_Q_match()/wilson_p("WPARAM_MATCH_SM", 6)/wilson_p("WPARAM_MATCH_SM", 6)))+950./243.+10./81.*PI*PI+124./27.*L+10./27.*L*L;
}

void C5::NNLO_calculation() {
    double L = wilson_p("WPARAM_MATCH_SM", 3);
    double coeff_temp = -G1t(wilson_p("WPARAM_MATCH_SM", {2,1}),log(this->get_Q_match()*this->get_Q_match()/wilson_p("WPARAM_MATCH_SM", 6)/wilson_p("WPARAM_MATCH_SM", 6)))/10.+2./15.*E0t(wilson_p("WPARAM_MATCH_SM", {2,1}))+68./243.+2./81.*PI*PI+14./81.*L+2./27.*L*L;
    // return this->double_to_complex_save("NNLO", coeff_temp);
}

void C6::NNLO_calculation() {
    double L = wilson_p("WPARAM_MATCH_SM", 3);
    double coeff_temp = -3./16.*G1t(wilson_p("WPARAM_MATCH_SM", {2,1}),log(this->get_Q_match()*this->get_Q_match()/wilson_p("WPARAM_MATCH_SM", 6)/wilson_p("WPARAM_MATCH_SM", 6)))+E0t(wilson_p("WPARAM_MATCH_SM", {2,1}))/4.+85./162.+5./108.*PI*PI+35./108.*L+5./36.*L*L;
    // return this->double_to_complex_save("NNLO", coeff_temp);
}

void C7::LO_calculation() {
    std::unordered_set<ParamId> sources {
        {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2, 1)}  // x_t
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        auto xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
        dep_param->set_expected(-0.5 * A0t(xt) - 23. / 36.);
    };

    WilsonParamComposer().compose_parameter(ParamId{ParameterType::WILSON, "B_MATCH", LhaID(305, 4422, 0, 0)}, sources, func);
}

void C7::NLO_calculation() {
    double L = wilson_p("WPARAM_MATCH_SM", 3);
    double coeff_temp = -0.5*A1t(wilson_p("WPARAM_MATCH_SM", {2,1}),log(this->get_Q_match()*this->get_Q_match()/wilson_p("WPARAM_MATCH_SM", 6)/wilson_p("WPARAM_MATCH_SM", 6)))+713./243.+4./81.*L-4./9.*(E0t(wilson_p("WPARAM_MATCH_SM", {2,1}))-7./9.+2./3.*L);;
    // return this->double_to_complex_save("NLO", coeff_temp);
}

void C7::NNLO_calculation() {
    double xtW=pow(QCDHelper::msbar_mass(6, sm("MASS",24))/sm("MASS", 24), 2); // mass top at pole for mtot param
	double xtt=pow(QCDHelper::mass_t_msbar()/sm("MASS",24),2.); // 24 -> W

    double coeff_temp = (C7t2mt(xtt)+log(this->get_Q_match()*this->get_Q_match()/wilson_p("WPARAM_MATCH_SM", 6)/wilson_p("WPARAM_MATCH_SM", 6))*((-592.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),5.)-22.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),4.)+12814.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),3.)-6376.*wilson_p("WPARAM_MATCH_SM", {2,1})*wilson_p("WPARAM_MATCH_SM", {2,1})+512.*wilson_p("WPARAM_MATCH_SM", {2,1}))/27./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,5.)*Li2(1.-1./wilson_p("WPARAM_MATCH_SM", {2,1}))
	+(-26838.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),5.)+25938.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),4.)+627367.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),3.)-331956.*wilson_p("WPARAM_MATCH_SM", {2,1})*wilson_p("WPARAM_MATCH_SM", {2,1})+16989.*wilson_p("WPARAM_MATCH_SM", {2,1})-460.)/729./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,6.)*log(wilson_p("WPARAM_MATCH_SM", {2,1}))
	+(34400.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),5.)+276644.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),4.)-2668324.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),3.)+1694437.*wilson_p("WPARAM_MATCH_SM", {2,1})*wilson_p("WPARAM_MATCH_SM", {2,1})-323354.*wilson_p("WPARAM_MATCH_SM", {2,1})+53077.)/2187./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,5.)
	+log(this->get_Q_match()*this->get_Q_match()/wilson_p("WPARAM_MATCH_SM", 6)/wilson_p("WPARAM_MATCH_SM", 6))*((-63.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),5.)+532.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),4.)+2089.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),3.)-1118.*wilson_p("WPARAM_MATCH_SM", {2,1})*wilson_p("WPARAM_MATCH_SM", {2,1}))/9./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,6.)*log(wilson_p("WPARAM_MATCH_SM", {2,1}))
	+(1186.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),5.)-2705.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),4.)-24791.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),3.)-16099.*wilson_p("WPARAM_MATCH_SM", {2,1})*wilson_p("WPARAM_MATCH_SM", {2,1})+19229.*wilson_p("WPARAM_MATCH_SM", {2,1})-2740.)/162./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,5.))) )
	-(C7c2MW(xtW)+13763./2187.*log(this->get_Q_match()*this->get_Q_match()/sm("MASS",24)/sm("MASS",24))+814./729.*pow(log(this->get_Q_match()*this->get_Q_match()/sm("MASS",24)/sm("MASS",24)),2.));
    // return this->double_to_complex_save("NNLO", coeff_temp);
}

void C8::LO_calculation() {
    std::unordered_set<ParamId> sources {
        {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2, 1)}  // x_t
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        auto xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
        dep_param->set_expected(-0.5*F0t(xt)-1./3.);
    };

    WilsonParamComposer().compose_parameter(ParamId{ParameterType::WILSON, "B_MATCH", LhaID(305, 6421, 0, 0)}, sources, func);


    // double coeff_temp = -0.5*F0t(wilson_p("WPARAM_MATCH_SM", {2,1}))-1./3.;
    // return this->double_to_complex_save("LO", coeff_temp);
}

void C8::NLO_calculation() {
    double L = wilson_p("WPARAM_MATCH_SM", 3);
    double coeff_temp = -0.5*F1t(wilson_p("WPARAM_MATCH_SM", {2,1}),log(this->get_Q_match()*this->get_Q_match()/wilson_p("WPARAM_MATCH_SM", 6)/wilson_p("WPARAM_MATCH_SM", 6)))+91./324.-4./27.*L-(E0t(wilson_p("WPARAM_MATCH_SM", {2,1}))-7./9.+2./3.*L)/6.;
    // return this->double_to_complex_save("NLO", coeff_temp);
}

void C8::NNLO_calculation() {

    double xtW=pow(QCDHelper::msbar_mass(6, sm("MASS",24))/sm("MASS", 24), 2); // mass top at pole for mtot param
	double xtt=pow(QCDHelper::mass_t_msbar()/sm("MASS",24),2.); // 24 -> W

    double coeff_temp = (C8t2mt(xtt)+log(this->get_Q_match()*this->get_Q_match()/wilson_p("WPARAM_MATCH_SM", 6)/wilson_p("WPARAM_MATCH_SM", 6))*((-148.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),5.)+1052.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),4.)-4811.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),3.)-3520.*wilson_p("WPARAM_MATCH_SM", {2,1})*wilson_p("WPARAM_MATCH_SM", {2,1})-61.*wilson_p("WPARAM_MATCH_SM", {2,1}))/18./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,5.)*Li2(1.-1./wilson_p("WPARAM_MATCH_SM", {2,1}))
	+(-15984.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),5.)+152379.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),4.)-1358060.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),3.)-1201653.*wilson_p("WPARAM_MATCH_SM", {2,1})*wilson_p("WPARAM_MATCH_SM", {2,1})-74190.*wilson_p("WPARAM_MATCH_SM", {2,1})+9188.)/1944./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,6.)*log(wilson_p("WPARAM_MATCH_SM", {2,1}))
	+(109669.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),5.)-1112675.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),4.)+6239377.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),3.)+8967623.*wilson_p("WPARAM_MATCH_SM", {2,1})*wilson_p("WPARAM_MATCH_SM", {2,1})+768722.*wilson_p("WPARAM_MATCH_SM", {2,1})-42796.)/11664./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,5.)
	+log(this->get_Q_match()*this->get_Q_match()/wilson_p("WPARAM_MATCH_SM", 6)/wilson_p("WPARAM_MATCH_SM", 6))*((-139.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),4.)-2938.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),3.)-2683.*wilson_p("WPARAM_MATCH_SM", {2,1})*wilson_p("WPARAM_MATCH_SM", {2,1}))/12./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,6.)*log(wilson_p("WPARAM_MATCH_SM", {2,1}))
	+(1295.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),5.)-7009.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),4.)+29495.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),3.)+64513.*wilson_p("WPARAM_MATCH_SM", {2,1})*wilson_p("WPARAM_MATCH_SM", {2,1})+17458.*wilson_p("WPARAM_MATCH_SM", {2,1})-2072.)/216./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,5.))) )
	-(C8c2MW(xtW)+16607./5832.*log(this->get_Q_match()*this->get_Q_match()/sm("MASS",24)/sm("MASS",24))+397./486.*pow(log(this->get_Q_match()*this->get_Q_match()/sm("MASS",24)/sm("MASS",24)),2.));
    // return this->double_to_complex_save("NNLO", coeff_temp);
}

void C9::LO_calculation() {
    std::unordered_set<ParamId> sources {
        {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2, 1)},  // x_t
        {ParameterType::WILSON, "WPARAM_SI_SM", 4},  // euh TODO
        {ParameterType::WILSON, "WPARAM_MATCH_SM", 3} //L
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double L = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 3})->get_val();
        auto xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
        double euh = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
        double coeff_temp = (1.-4.*euh)/euh*C0t(xt)-B0t(xt)/euh-D0t(xt) +1./4./euh+38./27.-4./9.*L;
        dep_param->set_expected(coeff_temp);
    };

    WilsonParamComposer().compose_parameter(ParamId{ParameterType::WILSON, "B_MATCH", LhaID(3051313, 4133, 0, 0)}, sources, func);


    // double L = wilson_p("WPARAM_MATCH_SM", 3);
    // double coeff_temp = (1.-4.*wilson_p("WPARAM_SI_SM", 4))/wilson_p("WPARAM_SI_SM", 4)*C0t(wilson_p("WPARAM_MATCH_SM", {2,1}))-B0t(wilson_p("WPARAM_MATCH_SM", {2,1}))/wilson_p("WPARAM_SI_SM", 4)-D0t(wilson_p("WPARAM_MATCH_SM", {2,1})) +1./4./wilson_p("WPARAM_SI_SM", 4)+38./27.-4./9.*L;
    // return this->double_to_complex_save("LO", coeff_temp);
}

void C9::NLO_calculation() {
    double L = wilson_p("WPARAM_MATCH_SM", 3);
    double coeff_temp = (1.-4.*wilson_p("WPARAM_SI_SM", 4))/wilson_p("WPARAM_SI_SM", 4)*C1t(wilson_p("WPARAM_MATCH_SM", {2,1}),log(this->get_Q_match()*this->get_Q_match()/wilson_p("WPARAM_MATCH_SM", 6)/wilson_p("WPARAM_MATCH_SM", 6)))
    -B1t(wilson_p("WPARAM_MATCH_SM", {2,1}),log(this->get_Q_match()*this->get_Q_match()/wilson_p("WPARAM_MATCH_SM", 6)/wilson_p("WPARAM_MATCH_SM", 6)))/wilson_p("WPARAM_SI_SM", 4)
    -D1t(wilson_p("WPARAM_MATCH_SM", {2,1}),log(this->get_Q_match()*this->get_Q_match()/wilson_p("WPARAM_MATCH_SM", 6)/wilson_p("WPARAM_MATCH_SM", 6))) +1./wilson_p("WPARAM_SI_SM", 4)+524./729.-128./243.*PI*PI-16./3.*L-128./81.*L*L;
    // return this->double_to_complex_save("NLO", coeff_temp);
}

void C9::NNLO_calculation() {
    this->set_CoefficientMatchingValue("LO", 1);
    // return 1;
}


void C10::LO_calculation() {

    std::unordered_set<ParamId> sources {
        {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2, 1)},  // x_t
        {ParameterType::WILSON, "WPARAM_SI_SM", 4}  // euh TODO
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        auto xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
        double euh = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
        double coeff_temp = (B0t(xt)-C0t(xt))/euh-1./4./euh;
        dep_param->set_expected(coeff_temp);
    };

    WilsonParamComposer().compose_parameter(ParamId{ParameterType::WILSON, "B_MATCH", LhaID(3051313, 4137, 0, 0)}, sources, func);


    double coeff_temp = (B0t(wilson_p("WPARAM_MATCH_SM", {2,1}))-C0t(wilson_p("WPARAM_MATCH_SM", {2,1})))/wilson_p("WPARAM_SI_SM", 4)-1./4./wilson_p("WPARAM_SI_SM", 4);
    // std::cout << coeff_temp << std::endl;
    // return this->double_to_complex_save("LO", coeff_temp);
}

void C10::NLO_calculation() {
    double coeff_temp =  (B1t(wilson_p("WPARAM_MATCH_SM", {2,1}),log(this->get_Q_match()*this->get_Q_match()/wilson_p("WPARAM_MATCH_SM", 6)/wilson_p("WPARAM_MATCH_SM", 6)))
    -C1t(wilson_p("WPARAM_MATCH_SM", {2,1}),log(this->get_Q_match()*this->get_Q_match()/wilson_p("WPARAM_MATCH_SM", 6)/wilson_p("WPARAM_MATCH_SM", 6))))/wilson_p("WPARAM_SI_SM", 4)-1./wilson_p("WPARAM_SI_SM", 4);
    // return this->double_to_complex_save("NLO", coeff_temp);
}
 
void C10::NNLO_calculation() {
    double xtW=pow(QCDHelper::msbar_mass(6, sm("MASS",24)) / sm("MASS", 24), 2); // mass top at pole for mtot param
	double xtt=pow(QCDHelper::mass_t_msbar()/sm("MASS",24),2.); // 24 -> W
    double coeff_temp = ((C10Wt2mt(xtt)+log(this->get_Q_match()*this->get_Q_match()/wilson_p("WPARAM_MATCH_SM", 6)/wilson_p("WPARAM_MATCH_SM", 6))*((69.+1292.*wilson_p("WPARAM_MATCH_SM", {2,1})-209.*wilson_p("WPARAM_MATCH_SM", {2,1})*wilson_p("WPARAM_MATCH_SM", {2,1}))/18./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,3.)
	-(521.*wilson_p("WPARAM_MATCH_SM", {2,1})+105.*wilson_p("WPARAM_MATCH_SM", {2,1})*wilson_p("WPARAM_MATCH_SM", {2,1})-50.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),3.))/9./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,4.)*log(wilson_p("WPARAM_MATCH_SM", {2,1}))
	-(47.*wilson_p("WPARAM_MATCH_SM", {2,1})+wilson_p("WPARAM_MATCH_SM", {2,1})*wilson_p("WPARAM_MATCH_SM", {2,1}))/3./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,3.)*Li2(1.-1./wilson_p("WPARAM_MATCH_SM", {2,1}))
	+log(this->get_Q_match()*this->get_Q_match()/wilson_p("WPARAM_MATCH_SM", 6)/wilson_p("WPARAM_MATCH_SM", 6))*((61.*wilson_p("WPARAM_MATCH_SM", {2,1})+11.*wilson_p("WPARAM_MATCH_SM", {2,1})*wilson_p("WPARAM_MATCH_SM", {2,1}))/3./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,3.)-(49.*wilson_p("WPARAM_MATCH_SM", {2,1})+96.*wilson_p("WPARAM_MATCH_SM", {2,1})*wilson_p("WPARAM_MATCH_SM", {2,1})-pow(wilson_p("WPARAM_MATCH_SM", {2,1}),3.))/6./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,4.)*log(wilson_p("WPARAM_MATCH_SM", {2,1})))))
	-(C10Wc2MW(xtW)-23./6.*log(this->get_Q_match()*this->get_Q_match()/sm("MASS", 24)/sm("MASS",24)))
	+(C10Zt2mt(xtt)+log(this->get_Q_match()*this->get_Q_match()/wilson_p("WPARAM_MATCH_SM", 6)/wilson_p("WPARAM_MATCH_SM", 6))*((188.*wilson_p("WPARAM_MATCH_SM", {2,1})+4.*wilson_p("WPARAM_MATCH_SM", {2,1})*wilson_p("WPARAM_MATCH_SM", {2,1})+95.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),3.)-47.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),4.))/6./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,3.)*Li2(1.-1./wilson_p("WPARAM_MATCH_SM", {2,1}))
	+(1468.*wilson_p("WPARAM_MATCH_SM", {2,1})+1578.*wilson_p("WPARAM_MATCH_SM", {2,1})*wilson_p("WPARAM_MATCH_SM", {2,1})-25.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),3.)-141.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),4.))/18./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,4.)*log(wilson_p("WPARAM_MATCH_SM", {2,1}))
	-(4622.*wilson_p("WPARAM_MATCH_SM", {2,1})+1031.*wilson_p("WPARAM_MATCH_SM", {2,1})*wilson_p("WPARAM_MATCH_SM", {2,1})+582.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),3.)-475.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),4.))/36./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,3.)
	+log(this->get_Q_match()*this->get_Q_match()/wilson_p("WPARAM_MATCH_SM", 6)/wilson_p("WPARAM_MATCH_SM", 6))*((49.*wilson_p("WPARAM_MATCH_SM", {2,1})+315.*wilson_p("WPARAM_MATCH_SM", {2,1})*wilson_p("WPARAM_MATCH_SM", {2,1})-4.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),3.))/6./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,4.)*log(wilson_p("WPARAM_MATCH_SM", {2,1}))-(440.*wilson_p("WPARAM_MATCH_SM", {2,1})+257.*wilson_p("WPARAM_MATCH_SM", {2,1})*wilson_p("WPARAM_MATCH_SM", {2,1})+72.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),3.)-49.*pow(wilson_p("WPARAM_MATCH_SM", {2,1}),4.))/12./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,3.))))
	+C10Z2tri(xtt)
	)*(-2./wilson_p("WPARAM_SI_SM", 4));
    // return this->double_to_complex_save("NNLO", coeff_temp);
}
void CQ1::LO_calculation() {
    double CSc_SM=-wilson_p("WPARAM_MATCH_SM", {2,1})*(wilson_p("WPARAM_MATCH_SM", {2,1})-2.)/12./(wilson_p("WPARAM_MATCH_SM", {2,1})-1.)/(wilson_p("WPARAM_MATCH_SM", {2,1})-1.)+(wilson_p("WPARAM_MATCH_SM", {2,1})-2.)*(3.*wilson_p("WPARAM_MATCH_SM", {2,1})-1.)/24./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,3.)*log(wilson_p("WPARAM_MATCH_SM", {2,1}));
    double CSn_SMonly=-3.*wilson_p("WPARAM_MATCH_SM", {2,1})/8./wilson_p("WPARAM_SI_SM", 1)+wilson_p("WPARAM_MATCH_SM", {2,1})*F0SP(wilson_p("WPARAM_MATCH_SM", {2,1}));
    double coeff_temp=(CSc_SM+CSn_SMonly)*(wilson_p("WPARAM_SI_SM", 3)*wilson_p("WPARAM_MATCH_SM", {5,2})/sm("MASS",24)/sm("MASS",24))/wilson_p("WPARAM_SI_SM", 4);
    // return this->double_to_complex_save("LO", coeff_temp);
}

void CQ2::LO_calculation() {
    double CPc_SM=1./24.*(wilson_p("WPARAM_MATCH_SM", {2,1})*(36.*wilson_p("WPARAM_MATCH_SM", {2,3})-203.*wilson_p("WPARAM_MATCH_SM", {2,2})+352.*wilson_p("WPARAM_MATCH_SM", {2,1})-209.)/6./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,3.)+(17.*wilson_p("WPARAM_MATCH_SM", {2,4})-34.*wilson_p("WPARAM_MATCH_SM", {2,3})+4.*wilson_p("WPARAM_MATCH_SM", {2,2})+23.*wilson_p("WPARAM_MATCH_SM", {2,1})-6.)/pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,4.)*log(wilson_p("WPARAM_MATCH_SM", {2,1})))
	-wilson_p("WPARAM_SI_SM", 4)/36.*(wilson_p("WPARAM_MATCH_SM", {2,1})*(18.*wilson_p("WPARAM_MATCH_SM", {2,3})-139.*wilson_p("WPARAM_MATCH_SM", {2,2})+274.*wilson_p("WPARAM_MATCH_SM", {2,1})-129.)/2./pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,3.)+(24.*wilson_p("WPARAM_MATCH_SM", {2,4})-33.*wilson_p("WPARAM_MATCH_SM", {2,3})-45.*wilson_p("WPARAM_MATCH_SM", {2,2})+50.*wilson_p("WPARAM_MATCH_SM", {2,1})-8.)/pow(wilson_p("WPARAM_MATCH_SM", {2,1})-1.,4.)*log(wilson_p("WPARAM_MATCH_SM", {2,1})));
    double CPn_SMonly=0.;
    double coeff_temp = (CPc_SM+CPn_SMonly)*(wilson_p("WPARAM_SI_SM", 3)*wilson_p("WPARAM_MATCH_SM", {5,2})/sm("MASS",24)/sm("MASS",24))/wilson_p("WPARAM_SI_SM", 4);
    // return this->double_to_complex_save("LO", coeff_temp);
}

void CP7::LO_calculation() {
    double coeff_temp = sm("MASS", 3) / wilson_p("WPARAM_MATCH_SM", {5,2}) * (-0.5 * A0t(wilson_p("WPARAM_MATCH_SM", {2,1})) - 23. / 36.);
    // return this->double_to_complex_save("LO", coeff_temp);
}

void CP8::LO_calculation() {
    double coeff_temp = sm("MASS", 3) / wilson_p("WPARAM_MATCH_SM", {5,2}) * (-0.5 * F0t(wilson_p("WPARAM_MATCH_SM", {2,1})) - 1. / 3.);
    // return this->double_to_complex_save("LO", coeff_temp);
}

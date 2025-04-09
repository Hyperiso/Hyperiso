#include "BWilson.h"

// ---------- C1 ----------

void C1::NLO_calculation() {
    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", 3}  // L
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        auto L = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 3})->get_val();
        dep_param->set_expected(15 + 6 * L);
    };

    WilsonParamComposer().compose_parameter(ParamId{"B_MATCH", LhaID(3040405, 6161, 1, 0)}, sources, func);
}

void C1::NNLO_calculation() {
    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", 3},              // L
        {"WPARAM_MATCH_SM", LhaID(2, 1)}     // x_t
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        auto L = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 3})->get_val();
        auto xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
        dep_param->set_expected(-T(xt) + 7987./72. + 17. * PI2 / 3. + 475./6. * L + 17. * L * L);
    };

    WilsonParamComposer().compose_parameter(ParamId{"B_MATCH", LhaID(3040405, 6161, 2, 0)}, sources, func);
}

// ---------- C2 ----------

void C2::LO_calculation() {
    std::unordered_set<ParamId> sources {};

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        dep_param->set_expected(1.);
    };

    WilsonParamComposer().compose_parameter(ParamId{"B_MATCH", LhaID(3040405, 4141, 0, 0)}, sources, func);
}

void C2::NNLO_calculation() {
    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", 3}  // L
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        auto L = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 3})->get_val();
        dep_param->set_expected(127. / 18. + 4. / 3. * PI2 + 46. / 3. * L + 4. * L * L);
    };

    WilsonParamComposer().compose_parameter(ParamId{"B_MATCH", LhaID(3040405, 4141, 2, 0)}, sources, func);
}

// ---------- C3 ----------

void C3::NNLO_calculation() {

    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", 3},  // L
        {"WPARAM_MATCH_SM", 6},
        {"WPARAM_MATCH_SM", {2, 1}},
        {"EW_SCALE", 1}
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double L = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 3})->get_val();
        double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
        double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();
        double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
        double coeff_temp = G1t(xt,log(Q_match*Q_match/mass_top_muW/mass_top_muW))-680./243.-20./81.*PI2-68./81.*L-20./27.*L*L;
        dep_param->set_expected(coeff_temp);
    };

    WilsonParamComposer().compose_parameter(ParamId{"B_MATCH", LhaID(3050707, 4133, 2, 0)}, sources, func);
}

// ---------- C4 ----------

void C4::NLO_calculation() {

    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", 3},  // L
        {"WPARAM_MATCH_SM", {2,1}}
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double L = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 3})->get_val();
        double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();
        double coeff_temp = E0t(xt)-7./9.+2./3.*L;
        dep_param->set_expected(coeff_temp);
    };

    WilsonParamComposer().compose_parameter(ParamId{"B_MATCH", LhaID(3050707, 6153, 1, 0)}, sources, func);
}

void C4::NNLO_calculation() {

    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", 3},  // L
        {"WPARAM_MATCH_SM", 6},
        {"WPARAM_MATCH_SM", {2,1}},
        {"EW_SCALE", 1}
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double L = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 3})->get_val();
        double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();
        double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
        double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
        double coeff_temp = E1t(xt,log(Q_match*Q_match/mass_top_muW/mass_top_muW))+950./243.+10./81.*PI2+124./27.*L+10./27.*L*L;
        dep_param->set_expected(coeff_temp);
    };

    WilsonParamComposer().compose_parameter(ParamId{"B_MATCH", LhaID(3050707, 6153, 2, 0)}, sources, func);
}

// ---------- C5 ----------

void C5::NNLO_calculation() {

    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", 3},  // L
        {"WPARAM_MATCH_SM", 6},
        {"WPARAM_MATCH_SM", {2,1}},
        {"EW_SCALE", 1}
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double L = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 3})->get_val();
        double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();
        double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
        double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
        double coeff_temp = -G1t(xt,log(Q_match*Q_match/mass_top_muW/mass_top_muW))/10.+2./15.*E0t(xt)+68./243.+2./81.*PI2+14./81.*L+2./27.*L*L;
        dep_param->set_expected(coeff_temp);
    };

    WilsonParamComposer().compose_parameter(ParamId{"B_MATCH", LhaID(3050707, 4536, 2, 0)}, sources, func);
}

// ---------- C6 ----------

void C6::NNLO_calculation() {

    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", 3},  // L
        {"WPARAM_MATCH_SM", 6},
        {"WPARAM_MATCH_SM", {2,1}},
        {"EW_SCALE", 1}
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double L = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 3})->get_val();
        double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();
        double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
        double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
        double coeff_temp = -3./16.*G1t(xt,log(Q_match*Q_match/mass_top_muW/mass_top_muW))+E0t(xt)/4.+85./162.+5./108.*PI2+35./108.*L+5./36.*L*L;
        dep_param->set_expected(coeff_temp);
    };

    WilsonParamComposer().compose_parameter(ParamId{"B_MATCH", LhaID(3050707, 6556, 2, 0)}, sources, func);
}

// ---------- C7 ----------

void C7::LO_calculation() {
    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", LhaID(2, 1)}  // x_t
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
        dep_param->set_expected(-0.5 * A0t(xt) - 23. / 36.);
    };

    WilsonParamComposer().compose_parameter(ParamId{"B_MATCH", LhaID(305, 4422, 0, 0)}, sources, func);
}

void C7::NLO_calculation() {

    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", 3},  // L
        {"WPARAM_MATCH_SM", 6},
        {"WPARAM_MATCH_SM", {2,1}},
        {"EW_SCALE", 1}
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double L = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 3})->get_val();
        double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();
        double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
        double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
        double coeff_temp = -0.5*A1t(xt,log(Q_match*Q_match/mass_top_muW/mass_top_muW))+713./243.+4./81.*L-4./9.*(E0t(xt)-7./9.+2./3.*L);
        dep_param->set_expected(coeff_temp);
    };

    WilsonParamComposer().compose_parameter(ParamId{"B_MATCH", LhaID(305, 4422, 1, 0)}, sources, func);
}

void C7::NNLO_calculation() {

    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", 3},  // L
        {"WPARAM_MATCH_SM", 6},
        {"WPARAM_MATCH_SM", {2,1}},
        {"WPARAM_MATCH_SM", 7},
        {"WPARAM_MATCH_SM", 8},
        {"EW_SCALE", 1},
        {ParameterType::SM, "MASS", 24}
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double L = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 3})->get_val();
        double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();
        double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
        double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();

        double xtW=src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 7})->get_val(); // mass top at pole for mtot param
	    double xtt=src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 8})->get_val(); // 24 -> W

        double m_W = src.at({ParameterType::SM, "MASS", 24})->get_val();
        double coeff_temp = (C7t2mt(xtt)+log(Q_match*Q_match/mass_top_muW/mass_top_muW)*((-592.*pow(xt,5.)-22.*pow(xt,4.)+12814.*pow(xt,3.)-6376.*xt*xt+512.*xt)/27./pow(xt-1.,5.)*Li2(1.-1./xt)
            +(-26838.*pow(xt,5.)+25938.*pow(xt,4.)+627367.*pow(xt,3.)-331956.*xt*xt+16989.*xt-460.)/729./pow(xt-1.,6.)*log(xt)
            +(34400.*pow(xt,5.)+276644.*pow(xt,4.)-2668324.*pow(xt,3.)+1694437.*xt*xt-323354.*xt+53077.)/2187./pow(xt-1.,5.)
            +log(Q_match*Q_match/mass_top_muW/mass_top_muW)*((-63.*pow(xt,5.)+532.*pow(xt,4.)+2089.*pow(xt,3.)-1118.*xt*xt)/9./pow(xt-1.,6.)*log(xt)
            +(1186.*pow(xt,5.)-2705.*pow(xt,4.)-24791.*pow(xt,3.)-16099.*xt*xt+19229.*xt-2740.)/162./pow(xt-1.,5.))) )
            -(C7c2MW(xtW)+13763./2187.*log(Q_match*Q_match/m_W/m_W)+814./729.*pow(log(Q_match*Q_match/m_W/m_W),2.));
            
            dep_param->set_expected(coeff_temp);
        };

    WilsonParamComposer().compose_parameter(ParamId{"B_MATCH", LhaID(305, 4422, 2, 0)}, sources, func);
}

// ---------- C8 ----------

void C8::LO_calculation() {
    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", LhaID(2, 1)}  // x_t
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
        dep_param->set_expected(-0.5*F0t(xt)-1./3.);
    };

    WilsonParamComposer().compose_parameter(ParamId{"B_MATCH", LhaID(305, 6421, 0, 0)}, sources, func);
}

void C8::NLO_calculation() {

    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", 3},  // L
        {"WPARAM_MATCH_SM", 6},
        {"WPARAM_MATCH_SM", {2,1}},
        {"EW_SCALE", 1},
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double L = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 3})->get_val();
        double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();
        double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
        double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();


        double coeff_temp = -0.5*F1t(xt,log(Q_match*Q_match/mass_top_muW/mass_top_muW))+91./324.-4./27.*L-(E0t(xt)-7./9.+2./3.*L)/6.;
        dep_param->set_expected(coeff_temp);
    };

    WilsonParamComposer().compose_parameter(ParamId{"B_MATCH", LhaID(305, 6421, 1, 0)}, sources, func);
}

void C8::NNLO_calculation() {

    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", 3},  // L
        {"WPARAM_MATCH_SM", 6},
        {"WPARAM_MATCH_SM", {2,1}},
        {"WPARAM_MATCH_SM", 7},
        {"WPARAM_MATCH_SM", 8},
        {"EW_SCALE", 1},
        {ParameterType::SM, "MASS", 24}
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double L = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 3})->get_val();
        double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();
        double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
        double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();

        double xtW=src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 7})->get_val(); // mass top at pole for mtot param
	    double xtt=src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 8})->get_val(); // 24 -> W

        double m_W = src.at({ParameterType::SM, "MASS", 24})->get_val();
        double coeff_temp = (C8t2mt(xtt)+log(Q_match*Q_match/mass_top_muW/mass_top_muW)*((-148.*pow(xt,5.)+1052.*pow(xt,4.)-4811.*pow(xt,3.)-3520.*xt*xt-61.*xt)/18./pow(xt-1.,5.)*Li2(1.-1./xt)
            +(-15984.*pow(xt,5.)+152379.*pow(xt,4.)-1358060.*pow(xt,3.)-1201653.*xt*xt-74190.*xt+9188.)/1944./pow(xt-1.,6.)*log(xt)
            +(109669.*pow(xt,5.)-1112675.*pow(xt,4.)+6239377.*pow(xt,3.)+8967623.*xt*xt+768722.*xt-42796.)/11664./pow(xt-1.,5.)
            +log(Q_match*Q_match/mass_top_muW/mass_top_muW)*((-139.*pow(xt,4.)-2938.*pow(xt,3.)-2683.*xt*xt)/12./pow(xt-1.,6.)*log(xt)
            +(1295.*pow(xt,5.)-7009.*pow(xt,4.)+29495.*pow(xt,3.)+64513.*xt*xt+17458.*xt-2072.)/216./pow(xt-1.,5.))) )
            -(C8c2MW(xtW)+16607./5832.*log(Q_match*Q_match/m_W/m_W)+397./486.*pow(log(Q_match*Q_match/m_W/m_W),2.));
        dep_param->set_expected(coeff_temp);
        };

    WilsonParamComposer().compose_parameter(ParamId{"B_MATCH", LhaID(305, 6421, 2, 0)}, sources, func);
}

// ---------- C9 ----------

void C9::LO_calculation() {
    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", LhaID(2, 1)},   // x_t
        {"WPARAM_MATCH_SM", 3},             // L
        {"WPARAM_SI_SM", 4},                // sw2
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double L = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 3})->get_val();
        double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
        double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
        double coeff_temp = (1.-4.*sw2)/sw2*C0t(xt)-B0t(xt)/sw2-D0t(xt) +1./4./sw2+38./27.-4./9.*L;
        dep_param->set_expected(coeff_temp);
    };

    WilsonParamComposer().compose_parameter(ParamId{"B_MATCH", LhaID(3051313, 4133, 0, 0)}, sources, func);
}

void C9::NLO_calculation() {

    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", 3},  // L
        {"WPARAM_MATCH_SM", 6},
        {"WPARAM_MATCH_SM", {2,1}},
        {"WPARAM_SI_SM", 4},
        {"EW_SCALE", 1}
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double L = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 3})->get_val();
        double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();
        double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
        double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();

        double sw2=src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val(); // sw2

        double coeff_temp = (1.-4.*sw2)/sw2*C1t(xt,log(Q_match*Q_match/mass_top_muW/mass_top_muW))
            -B1t(xt,log(Q_match*Q_match/mass_top_muW/mass_top_muW))/sw2
            -D1t(xt,log(Q_match*Q_match/mass_top_muW/mass_top_muW)) +1./sw2+524./729.-128./243.*PI2-16./3.*L-128./81.*L*L;
        dep_param->set_expected(coeff_temp);
        };

    WilsonParamComposer().compose_parameter(ParamId{"B_MATCH", LhaID(3051313, 4133, 1, 0)}, sources, func);
}

void C9::NNLO_calculation() {
    // this->set_WilsonCoeffMatching("LO", 1);
    // return 1;
}

// ---------- C10 ----------

void C10::LO_calculation() {
    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", LhaID(2, 1)},   // x_t
        {"WPARAM_SI_SM", 4}                 // sw2
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
        double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
        double coeff_temp = (B0t(xt) - C0t(xt) - .25) / sw2;
        dep_param->set_expected(coeff_temp);
    };

    WilsonParamComposer().compose_parameter(ParamId{"B_MATCH", LhaID(3051313, 4137, 0, 0)}, sources, func);
}

void C10::NLO_calculation() {

    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", 6},
        {"WPARAM_MATCH_SM", {2,1}},
        {"WPARAM_SI_SM", 4},
        {"EW_SCALE", 1}
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();
        double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
        double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();

        double sw2=src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val(); // sw2

        double coeff_temp =  (B1t(xt,log(Q_match*Q_match/mass_top_muW/mass_top_muW))
            -C1t(xt,log(Q_match*Q_match/mass_top_muW/mass_top_muW)))/sw2-1./sw2;
        dep_param->set_expected(coeff_temp);
        };

    WilsonParamComposer().compose_parameter(ParamId{"B_MATCH", LhaID(3051313, 4137, 1, 0)}, sources, func);
}
 
void C10::NNLO_calculation() {

    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", 6},
        {"WPARAM_MATCH_SM", {2,1}},
        {"WPARAM_MATCH_SM", 7},
        {"WPARAM_MATCH_SM", 8},
        {"WPARAM_SI_SM", 4},
        {"EW_SCALE", 1},
        {ParameterType::SM, "MASS", 24}
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();
        double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
        double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();

        double xtW=src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 7})->get_val(); // mass top at pole for mtot param
	    double xtt=src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 8})->get_val(); // 24 -> W

        double sw2=src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val(); // sw2

        double m_W = src.at({ParameterType::SM, "MASS", 24})->get_val();
        double coeff_temp = ((C10Wt2mt(xtt)+log(Q_match*Q_match/mass_top_muW/mass_top_muW)*((69.+1292.*xt-209.*xt*xt)/18./pow(xt-1.,3.)
            -(521.*xt+105.*xt*xt-50.*pow(xt,3.))/9./pow(xt-1.,4.)*log(xt)
            -(47.*xt+xt*xt)/3./pow(xt-1.,3.)*Li2(1.-1./xt)
            +log(Q_match*Q_match/mass_top_muW/mass_top_muW)*((61.*xt+11.*xt*xt)/3./pow(xt-1.,3.)-(49.*xt+96.*xt*xt-pow(xt,3.))/6./pow(xt-1.,4.)*log(xt))))
            -(C10Wc2MW(xtW)-23./6.*log(Q_match*Q_match/m_W/m_W))
            +(C10Zt2mt(xtt)+log(Q_match*Q_match/mass_top_muW/mass_top_muW)*((188.*xt+4.*xt*xt+95.*pow(xt,3.)-47.*pow(xt,4.))/6./pow(xt-1.,3.)*Li2(1.-1./xt)
            +(1468.*xt+1578.*xt*xt-25.*pow(xt,3.)-141.*pow(xt,4.))/18./pow(xt-1.,4.)*log(xt)
            -(4622.*xt+1031.*xt*xt+582.*pow(xt,3.)-475.*pow(xt,4.))/36./pow(xt-1.,3.)
            +log(Q_match*Q_match/mass_top_muW/mass_top_muW)*((49.*xt+315.*xt*xt-4.*pow(xt,3.))/6./pow(xt-1.,4.)*log(xt)-(440.*xt+257.*xt*xt+72.*pow(xt,3.)-49.*pow(xt,4.))/12./pow(xt-1.,3.))))
            +C10Z2tri(xtt)
            )*(-2./sw2);
        dep_param->set_expected(coeff_temp);
        };

    WilsonParamComposer().compose_parameter(ParamId{"B_MATCH", LhaID(3051313, 4137, 2, 0)}, sources, func);
}

// ---------- CQ1 ----------

void CQ1::LO_calculation() {

    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", {2,1}},
        {"WPARAM_SI_SM", 4},
        {"WPARAM_SI_SM", 3},
        {"WPARAM_SI_SM", 1},
        {ParameterType::SM, "MASS", 24}
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();
        double mass_b_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5,2}})->get_val(); // mass top at pole for mtot param

        double xh=src.at({ParameterType::WILSON, "WPARAM_SI_SM", 1})->get_val(); // xh
        double ml=src.at({ParameterType::WILSON, "WPARAM_SI_SM", 3})->get_val(); // ml
        double sw2=src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val(); // sw2

        double m_W = src.at({ParameterType::SM, "MASS", 24})->get_val();
        double CSc_SM=-xt*(xt-2.)/12./(xt-1.)/(xt-1.)+(xt-2.)*(3.*xt-1.)/24./pow(xt-1.,3.)*log(xt);
        double CSn_SMonly=-3.*xt/8./xh+xt*F0SP(xt);
        double coeff_temp=(CSc_SM+CSn_SMonly)*(ml*mass_b_muW/m_W/m_W)/sw2;
        dep_param->set_expected(coeff_temp);
    };

    WilsonParamComposer().compose_parameter(ParamId{"B_SCALAR_MATCH", LhaID(3051313, 3230, 0, 0)}, sources, func);
}

// ---------- CQ2 ----------

void CQ2::LO_calculation() {

    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", {2,1}},
        {"WPARAM_MATCH_SM", {2,2}},
        {"WPARAM_MATCH_SM", {2,3}},
        {"WPARAM_MATCH_SM", {2,4}},
        {"WPARAM_SI_SM", 4},
        {"WPARAM_SI_SM", 3},
        {ParameterType::SM, "MASS", 24}
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();
        double xt2 = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,2}})->get_val();
        double xt3 = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,3}})->get_val();
        double xt4 = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,4}})->get_val();
        double mass_b_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}})->get_val();

        double ml=src.at({ParameterType::WILSON, "WPARAM_SI_SM", 3})->get_val(); // ml
        double sw2=src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val(); // sw2

        double m_W = src.at({ParameterType::SM, "MASS", 24})->get_val();
        double CPc_SM=1./24.*(xt*(36.*xt3-203.*xt2+352.*xt-209.)/6./pow(xt-1.,3.)+(17.*xt4-34.*xt3+4.*xt2+23.*xt-6.)/pow(xt-1.,4.)*log(xt))
            -sw2/36.*(xt*(18.*xt3-139.*xt2+274.*xt-129.)/2./pow(xt-1.,3.)+(24.*xt4-33.*xt3-45.*xt2+50.*xt-8.)/pow(xt-1.,4.)*log(xt));
        double CPn_SMonly=0.;
        double coeff_temp = (CPc_SM+CPn_SMonly)*(ml*mass_b_muW/m_W/m_W)/sw2;
        dep_param->set_expected(coeff_temp);
    };

    WilsonParamComposer().compose_parameter(ParamId{"B_SCALAR_MATCH", LhaID(3051313, 3233, 0, 0)}, sources, func);
}

// ---------- C'7 ----------

void CP7::LO_calculation() {

    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", {2,1}},
        {"WPARAM_MATCH_SM", {5,1}},
        {"WPARAM_SI_SM", 4},
        {ParameterType::SM, "MASS", 3}
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();
        double mass_b_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}})->get_val();
        double m_c = src.at({ParameterType::SM, "MASS", 3})->get_val();
        double coeff_temp = m_c / mass_b_muW * (-0.5 * A0t(xt) - 23. / 36.);
        dep_param->set_expected(coeff_temp);
        };

    WilsonParamComposer().compose_parameter(ParamId{"B_PRIME_MATCH", LhaID(305, 4322, 0, 0)}, sources, func);
}

// ---------- C'8 ----------

void CP8::LO_calculation() {

    std::unordered_set<ParamId> sources {
        {"WPARAM_MATCH_SM", {2,1}},
        {"WPARAM_MATCH_SM", {5,1}},
        {"WPARAM_SI_SM", 4},
        {ParameterType::SM, "MASS", 3}
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();
        double mass_b_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}})->get_val();

        double m_c = src.at({ParameterType::SM, "MASS", 3})->get_val();
        double coeff_temp = m_c / mass_b_muW * (-0.5 * F0t(xt) - 1. / 3.);
        dep_param->set_expected(coeff_temp);
        };

    WilsonParamComposer().compose_parameter(ParamId{"B_PRIME_MATCH", LhaID(305, 4321, 0, 0)}, sources, func);
}

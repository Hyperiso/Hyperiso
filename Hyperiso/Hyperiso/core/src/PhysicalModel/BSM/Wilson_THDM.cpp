// #include "Wilson_THDM.h"

// //TODO ALL FUNCTION AND RETURN

// void C3_THDM::NNLO_calculation() {
//     std::unordered_set<ParamId> sources {
//         {"WPARAM_SI_BSM", 7},
//         {"WPARAM_MATCH_BSM", 1},
//         {"EW_SCALE", 1},
//         {ParameterType::BSM, "MASS", 37}
//     };

//     auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
//         auto yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
//         double m_H = src.at({ParameterType::BSM, "MASS", 37})->get_val();
//         double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
//         double coeff_temp = G3H(yt,lu)+Delta3H(yt,lu)*log(pow(Q_match/m_H,2.));
//         dep_param->set_expected(coeff_temp);
//     };

//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3050707, 4133, 2, 1)}, sources, func);

//     // double coeff_temp = G3H(thdm_params->yt,thdm_params->lu)+Delta3H(thdm_params->yt,thdm_params->lu)*log(pow(this->get_Q_match()/(*mod)("MASS",37),2.));
//     // // return this->double_to_complex_save("NNLO", coeff_temp);
// }

// void C4_THDM::NLO_calculation() {
//     std::unordered_set<ParamId> sources {
//         {"WPARAM_SI_BSM", 7},
//         {"WPARAM_MATCH_BSM", 1},
//     };

//     auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
//         auto yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
//         double coeff_temp = EH(yt,lu);
//         dep_param->set_expected(coeff_temp);
//     };

//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3050707, 6153, 1, 1)}, sources, func);

//     // double coeff_temp = EH(thdm_params->yt,thdm_params->lu);
//     // // return this->double_to_complex_save("NLO", coeff_temp);
// }

// void C4_THDM::NNLO_calculation() {
//     std::unordered_set<ParamId> sources {
//         {"WPARAM_SI_BSM", 7},
//         {"WPARAM_MATCH_BSM", 1},
//         {"EW_SCALE", 1},
//         {ParameterType::BSM, "MASS", 37}
//     };

//     auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
//         auto yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
//         double m_H = src.at({ParameterType::BSM, "MASS", 37})->get_val();
//         double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();

//         double coeff_temp=G4H(yt,lu)+Delta4H(yt,lu)*log(pow(Q_match/m_H,2.));
//         dep_param->set_expected(coeff_temp);
//     };

//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3050707, 6153, 2, 1)}, sources, func);

//     // double coeff_temp=G4H(thdm_params->yt,thdm_params->lu)+Delta4H(thdm_params->yt,thdm_params->lu)*log(pow(this->get_Q_match()/(*mod)("MASS",37),2.));
//     // // return this->double_to_complex_save("NNLO", coeff_temp);
// }

// void C5_THDM::NNLO_calculation() {

//     std::unordered_set<ParamId> sources {
//         {"WPARAM_SI_BSM", 7},
//         {"WPARAM_MATCH_BSM", 1},
//         {"EW_SCALE", 1},
//         {ParameterType::BSM, "MASS", 37}
//     };

//     auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
//         auto yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
//         double m_H = src.at({ParameterType::BSM, "MASS", 37})->get_val();
//         double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();

//         double C4H_1=EH(yt,lu);
//         double C3H_2=G3H(yt,lu)+Delta3H(yt,lu)*log(pow(Q_match/m_H,2.));
//         double coeff_temp=-C3H_2/10.+2./15.*C4H_1;
//         dep_param->set_expected(coeff_temp);
//     };

//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3050707, 4536, 2, 1)}, sources, func);


//     // double C4H_1=EH(thdm_params->yt,thdm_params->lu);
//     // double C3H_2=G3H(thdm_params->yt,thdm_params->lu)+Delta3H(thdm_params->yt,thdm_params->lu)*log(pow(this->get_Q_match()/(*mod)("MASS",37),2.));
// 	// double coeff_temp=-C3H_2/10.+2./15.*C4H_1;
//     // // return this->double_to_complex_save("NNLO", coeff_temp);
// }

// void C6_THDM::NNLO_calculation() {

//     std::unordered_set<ParamId> sources {
//         {"WPARAM_SI_BSM", 7},
//         {"WPARAM_MATCH_BSM", 1},
//         {"EW_SCALE", 1},
//         {ParameterType::BSM, "MASS", 37}
//     };

//     auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
//         auto yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
//         double m_H = src.at({ParameterType::BSM, "MASS", 37})->get_val();
//         double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();

//         double C4H_1=EH(yt,lu);
//         double C3H_2=G3H(yt,lu)+Delta3H(yt,lu)*log(pow(Q_match/m_H,2.));
//         double coeff_temp=-3./16.*C3H_2+1./4.*C4H_1;
//         dep_param->set_expected(coeff_temp);
//     };

//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3050707, 6556, 2, 1)}, sources, func);


//     // double C4H_1=EH(thdm_params->yt,thdm_params->lu);
//     // double C3H_2=G3H(thdm_params->yt,thdm_params->lu)+Delta3H(thdm_params->yt,thdm_params->lu)*log(pow(this->get_Q_match()/(*mod)("MASS",37),2.));
// 	// double coeff_temp=-3./16.*C3H_2+1./4.*C4H_1;

//     // // return this->double_to_complex_save("NNLO", coeff_temp);
// }

// void C7_THDM::LO_calculation() {

//     std::unordered_set<ParamId> sources {
//         {"WPARAM_SI_BSM", 7},
//         {"WPARAM_SI_BSM", 8},
//         {"WPARAM_MATCH_BSM", 1}
//     };

//     auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
//         double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
//         auto yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();

//         double coeff_temp = 1./3.*lu*lu*F7_1(yt) - lu*ld*F7_2(yt);
//         dep_param->set_expected(coeff_temp);
//     };

//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(305, 4422, 0, 1)}, sources, func);

//     // std::cout << "YUUU" <<  (*mod)("YU", 22) << std::endl;
//     // double coeff_temp = 1./3.*thdm_params->lu*thdm_params->lu*F7_1(thdm_params->yt) - thdm_params->lu*thdm_params->ld*F7_2(thdm_params->yt);
//     // // std::cout << "lu : " << thdm_params->lu << std::endl;
//     // // std::cout << "yt : " << thdm_params->yt << std::endl;
//     // // std::cout << "ld : " << thdm_params->ld << std::endl;
//     // // return this->double_to_complex_save("LO", coeff_temp);
// }

// void C7_THDM::NLO_calculation() {

//     std::unordered_set<ParamId> sources {
//         {"WPARAM_SI_BSM", 7},
//         {"WPARAM_SI_BSM", 8},
//         {"WPARAM_MATCH_BSM", 1},
//         {"EW_SCALE", 1},
//         {ParameterType::BSM, "MASS", 37}
//     };

//     auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
//         double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
//         auto yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
//         double m_H = src.at({ParameterType::BSM, "MASS", 37})->get_val();
//         double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();

//         double coeff_temp = G7H(yt,lu,ld)+Delta7H(yt,lu,ld)*log(pow(Q_match/m_H,2.))
//             -4./9.*EH(yt,lu);
//         dep_param->set_expected(coeff_temp);
//     };

//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(305, 4422, 1, 1)}, sources, func);


//     // double coeff_temp = G7H(thdm_params->yt,thdm_params->lu,thdm_params->ld)+Delta7H(thdm_params->yt,thdm_params->lu,thdm_params->ld)*log(pow(this->get_Q_match()/thdm_params->m_H,2.))
//     // -4./9.*EH(thdm_params->yt,thdm_params->lu);
//     // // return this->double_to_complex_save("NLO", coeff_temp);
// }

// void C7_THDM::NNLO_calculation() {

//     std::unordered_set<ParamId> sources {
//         {"WPARAM_SI_BSM", 7},
//         {"WPARAM_SI_BSM", 8},
//         {"WPARAM_MATCH_BSM", 1},
//         {"WPARAM_MATCH_SM", 6},
//         {"EW_SCALE", 1}
//     };

//     auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
//         double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
//         auto yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
//         double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
//         double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();

//         double coeff_temp =C7H2(yt,lu,ld,log(pow(Q_match/mass_top_muW, 2.)));
//         dep_param->set_expected(coeff_temp);
//     };

//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(305, 4422, 2, 1)}, sources, func);

// 	// double coeff_temp =C7H2(thdm_params->yt,thdm_params->lu,thdm_params->ld,log(pow(this->get_Q_match()/thdm_params->mass_top_muW, 2.)));
//     // // return this->double_to_complex_save("NNLO", coeff_temp);
// }

// void C8_THDM::LO_calculation() {

//     std::unordered_set<ParamId> sources {
//         {"WPARAM_SI_BSM", 7},
//         {"WPARAM_SI_BSM", 8},
//         {"WPARAM_MATCH_BSM", 1}
//     };

//     auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
//         double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
//         auto yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();

//         double coeff_temp = 1./3.*lu*ld*F8_1(yt) - lu*ld*F8_2(yt);
//         dep_param->set_expected(coeff_temp);
//     };

//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(305, 6421, 0, 1)}, sources, func);

//     // double coeff_temp = 1./3.*thdm_params->lu*thdm_params->lu*F8_1(thdm_params->yt) - thdm_params->lu*thdm_params->ld*F8_2(thdm_params->yt);
//     // // return this->double_to_complex_save("LO", coeff_temp);
// }

// void C8_THDM::NLO_calculation() {

//     std::unordered_set<ParamId> sources {
//         {"WPARAM_SI_BSM", 7},
//         {"WPARAM_SI_BSM", 8},
//         {"WPARAM_MATCH_BSM", 1},
//         {"EW_SCALE", 1},
//         {ParameterType::BSM, "MASS", 37}
//     };

//     auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
//         double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
//         auto yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
//         double m_H = src.at({ParameterType::BSM, "MASS", 37})->get_val();
//         double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();

//         double coeff_temp = G8H(yt,lu,ld)+Delta8H(yt,lu,ld)*log(pow(Q_match/m_H,2.))
//             -1./6.*EH(yt,lu);
//         dep_param->set_expected(coeff_temp);
//     };

//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(305, 6421, 1, 1)}, sources, func);


//     // double coeff_temp = G8H(thdm_params->yt,thdm_params->lu,thdm_params->ld)+Delta8H(thdm_params->yt,thdm_params->lu,thdm_params->ld)*log(pow(this->get_Q_match()/thdm_params->m_H,2.))
//     // -1./6.*EH(thdm_params->yt,thdm_params->lu);
//     // // return this->double_to_complex_save("NLO", coeff_temp);
// }

// void C8_THDM::NNLO_calculation() {

//     std::unordered_set<ParamId> sources {
//         {"WPARAM_SI_BSM", 7},
//         {"WPARAM_SI_BSM", 8},
//         {"WPARAM_MATCH_BSM", 1},
//         {"WPARAM_MATCH_SM", 6},
//         {"EW_SCALE", 1}
//     };

//     auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
//         double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
//         auto yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
//         double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
//         double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();

//         double coeff_temp =C8H2(yt,lu,ld,log(pow(Q_match/mass_top_muW, 2.)));
//         dep_param->set_expected(coeff_temp);
//     };

//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(305, 6421, 1, 1)}, sources, func);


// 	// double coeff_temp =C8H2(thdm_params->yt,thdm_params->lu,thdm_params->ld,log(pow(this->get_Q_match()/thdm_params->mass_top_muW, 2.)));
//     // // return this->double_to_complex_save("NNLO", coeff_temp);
// }

// void C9_THDM::LO_calculation() {

//     std::unordered_set<ParamId> sources {
//         {"WPARAM_SI_SM", 4},
//         {"WPARAM_SI_BSM", 7},
//         {"WPARAM_MATCH_BSM", 1},
//         {"WPARAM_MATCH_SM", {2,1}}
//     };

//     auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
//         double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
//         auto yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
//         auto xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();

//         double coeff_temp = (1.-4.*sw2)/sw2*C9llH0(xt,yt,lu)-D9H0(yt,lu);
//         dep_param->set_expected(coeff_temp);
//     };

//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3051313, 4133, 0, 1)}, sources, func);


//     // double coeff_temp = (1.-4.*thdm_params->sw2)/thdm_params->sw2*C9llH0(thdm_params->xt,thdm_params->yt,thdm_params->lu)-D9H0(thdm_params->yt,thdm_params->lu);
//     // // return this->double_to_complex_save("LO", coeff_temp);
// }

// void C9_THDM::NLO_calculation() {

//     std::unordered_set<ParamId> sources {
//         {"WPARAM_SI_SM", 4},
//         {"WPARAM_SI_BSM", 7},
//         {"WPARAM_MATCH_BSM", 1},
//         {"WPARAM_MATCH_SM", {2,1}},
//         {"EW_SCALE", 1},
//         {ParameterType::BSM, "MASS", 37}
//     };

//     auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
//         double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
//         auto yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
//         auto xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();
//         double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
//         double m_H = src.at({ParameterType::BSM, "MASS", 37})->get_val();

//         double coeff_temp = (1.-4.*sw2)/sw2*C9llH1(xt,yt,lu,log(pow(Q_match/m_H,2.)))
//             -D9H1(yt,lu,log(pow(Q_match/m_H,2.)));
//         dep_param->set_expected(coeff_temp);
//     };

//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3051313, 4133, 1, 1)}, sources, func);


//     // double coeff_temp = (1.-4.*thdm_params->sw2)/thdm_params->sw2*C9llH1(thdm_params->xt,thdm_params->yt,thdm_params->lu,log(pow(this->get_Q_match()/thdm_params->m_H,2.)))
//     // -D9H1(thdm_params->yt,thdm_params->lu,log(pow(this->get_Q_match()/thdm_params->m_H,2.)));
//     // // return this->double_to_complex_save("NLO", coeff_temp);
// }

// void C10_THDM::LO_calculation() {

//     std::unordered_set<ParamId> sources {
//         {"WPARAM_SI_SM", 4},
//         {"WPARAM_SI_BSM", 7},
//         {"WPARAM_MATCH_BSM", 1},
//         {"WPARAM_MATCH_SM", {2,1}}
//     };

//     auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
//         double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
//         auto yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
//         auto xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();

//         double coeff_temp = -C9llH0(xt,yt,lu)/sw2;
//         dep_param->set_expected(coeff_temp);
//     };

//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3051313, 4137, 0, 1)}, sources, func);

//     // double coeff_temp = -C9llH0(thdm_params->xt,thdm_params->yt,thdm_params->lu)/thdm_params->sw2;
//     // // return this->double_to_complex_save("LO", coeff_temp);
// }

// void C10_THDM::NLO_calculation() {

//     std::unordered_set<ParamId> sources {
//         {"WPARAM_SI_SM", 4},
//         {"WPARAM_SI_BSM", 7},
//         {"WPARAM_MATCH_BSM", 1},
//         {"WPARAM_MATCH_SM", {2,1}},
//         {"EW_SCALE", 1},
//         {ParameterType::BSM, "MASS", 37}
//     };

//     auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
//         double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
//         auto yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
//         auto xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();
//         double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
//         double m_H = src.at({ParameterType::BSM, "MASS", 37})->get_val();

//         double coeff_temp = -C9llH1(xt,yt,lu,log(pow(Q_match/m_H,2.)))/sw2;
//         dep_param->set_expected(coeff_temp);
//     };

//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3051313, 4137, 1, 1)}, sources, func);

//     // double coeff_temp = -C9llH1(thdm_params->xt,thdm_params->yt,thdm_params->lu,log(pow(this->get_Q_match()/thdm_params->m_H,2.)))/thdm_params->sw2;
//     // // return this->double_to_complex_save("NLO", coeff_temp);
// }

// void CQ1_THDM::LO_calculation() {

//     std::unordered_set<ParamId> sources {
//         {"WPARAM_SI_SM", 3},
//         {"WPARAM_SI_SM", 4},
//         {"WPARAM_SI_BSM", 2},
//         {"WPARAM_SI_BSM", 3},
//         {"WPARAM_SI_BSM", 6},
//         {"WPARAM_SI_BSM", 7},
//         {"WPARAM_SI_BSM", 8},
//         {"WPARAM_SI_BSM", 9},
//         {"WPARAM_SI_BSM", 10},
//         {"WPARAM_MATCH_SM", {2,1}},
//         {"WPARAM_MATCH_SM", 6},
//         {ParameterType::SM, "MASS", 24}
//     };

//     auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
//         double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
//         double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
        
//         double alpha = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 9})->get_val();
//         double beta = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val();

//         double xH = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 2})->get_val();
//         double xH0 = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 3})->get_val();
//         double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();
//         double m_W = src.at({ParameterType::SM, "MASS", 24})->get_val();
//         double ml = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 3})->get_val();
//         double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
//         double le = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 10})->get_val();
//         double G1=-3./4.+ld*lu*F4SP(xt,xH)+lu*lu*F5SP(xt,xH);
//         double G2=ld*(ld*lu+1.)*F6SP(xt,xH)-ld*lu*lu*F7SP(xt,xH)
//             +lu*lu*(ld*F8SP(xt,xH)+lu*F9SP(xt,xH)-lu*F10SP(xt,xH))+lu*F11SP(xt,xH)-lu*F12SP(xt,xH);

//         double CSn_2HDM=xt*(F0SP(xt)+le*(ld*F1SP(xt,xH)+lu*F2SP(xt,xH))+le*lu*F3SP(xt,xH))
//             +xt/2./xH*(sin(alpha-beta)+cos(alpha-beta)*le)*(sin(alpha-beta)*G1+cos(alpha-beta)*G2)
//             +xt/2./xH0*(cos(alpha-beta)-sin(alpha-beta)*le)*(cos(alpha-beta)*G1-sin(alpha-beta)*G2);
//         double coeff_temp=CSc_2HDM(xH,xt,lu,ld,le)+CSn_2HDM;
//         coeff_temp*=(ml*mass_top_muW/m_W/m_W)/sw2;
//         dep_param->set_expected(coeff_temp);
//     };

//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3051313, 3230, 0, 1)}, sources, func);

//     // double le = (*mod)("YL",10*(gen-1)+gen-1);
// 	// double G1=-3./4.+thdm_params->ld*thdm_params->lu*F4SP(thdm_params->xt,thdm_params->xH)+thdm_params->lu*thdm_params->lu*F5SP(thdm_params->xt,thdm_params->xH);
// 	// double G2=thdm_params->ld*(thdm_params->ld*thdm_params->lu+1.)*F6SP(thdm_params->xt,thdm_params->xH)-thdm_params->ld*thdm_params->lu*thdm_params->lu*F7SP(thdm_params->xt,thdm_params->xH)
// 	// +thdm_params->lu*thdm_params->lu*(thdm_params->ld*F8SP(thdm_params->xt,thdm_params->xH)+thdm_params->lu*F9SP(thdm_params->xt,thdm_params->xH)-thdm_params->lu*F10SP(thdm_params->xt,thdm_params->xH))+thdm_params->lu*F11SP(thdm_params->xt,thdm_params->xH)-thdm_params->lu*F12SP(thdm_params->xt,thdm_params->xH);

// 	// double CSn_2HDM=thdm_params->xt*(F0SP(thdm_params->xt)+le*(thdm_params->ld*F1SP(thdm_params->xt,thdm_params->xH)+thdm_params->lu*F2SP(thdm_params->xt,thdm_params->xH))+le*thdm_params->lu*F3SP(thdm_params->xt,thdm_params->xH))
// 	// +thdm_params->xt/2./thdm_params->xh*(sin(thdm_params->alpha-thdm_params->beta)+cos(thdm_params->alpha-thdm_params->beta)*le)*(sin(thdm_params->alpha-thdm_params->beta)*G1+cos(thdm_params->alpha-thdm_params->beta)*G2)
// 	// +thdm_params->xt/2./thdm_params->xH0*(cos(thdm_params->alpha-thdm_params->beta)-sin(thdm_params->alpha-thdm_params->beta)*le)*(cos(thdm_params->alpha-thdm_params->beta)*G1-sin(thdm_params->alpha-thdm_params->beta)*G2);
// 	// double coeff_temp=CSc_2HDM(thdm_params->xH,thdm_params->xt,thdm_params->lu,thdm_params->ld,le)+CSn_2HDM;
// 	// coeff_temp*=(wilson_p("WPARAM_SI_SM", 3)*thdm_params->mass_b_muW/sm("MASS",24)/sm("MASS",24))/thdm_params->sw2;
// 	// // return this->double_to_complex_save("LO", coeff_temp);
// }

// void CQ2_THDM::LO_calculation() {

//     std::unordered_set<ParamId> sources {
//         {"WPARAM_SI_SM", 3},
//         {"WPARAM_SI_SM", 4},
//         {"WPARAM_SI_BSM", 2},
//         {"WPARAM_SI_BSM", 4},
//         {"WPARAM_SI_BSM", 7},
//         {"WPARAM_SI_BSM", 8},
//         {"WPARAM_SI_BSM", 10},
//         {"WPARAM_MATCH_SM", {2,1}},
//         {"WPARAM_MATCH_SM", {5,1}},
//         {ParameterType::SM, "MASS", 24}
//     };

//     auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
//         double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
//         double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();

//         double xH = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 2})->get_val();
//         double xA = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 4})->get_val();
//         double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();
//         double m_W = src.at({ParameterType::SM, "MASS", 24})->get_val();
//         double ml = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 3})->get_val();
//         double mass_b_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}})->get_val();
//         double le = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 10})->get_val();

//         double G3=ld*(ld*lu+1.)*F6SP(xt,xH)+ld*lu*lu*F7SP(xt,xH)
//             +lu*lu*(ld*F8SP(xt,xH)+lu*F9SP(xt,xH)+lu*F10SP(xt,xH))+lu*F11SP(xt,xH)+lu*F12SP(xt,xH);
//         double CPn_2HDM=xt*(-le*(ld*F1SP(xt,xH)+lu*F2SP(xt,xH))+le*lu*F3SP(xt,xH))+xt/2./xA*(le)*G3;

//         double coeff_temp=CPc_2HDM(xH,xt,lu,ld,le,sw2)+CPn_2HDM;
//         coeff_temp*=(ml*mass_b_muW/m_W/m_W)/sw2;
//         dep_param->set_expected(coeff_temp);
//     };

//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3051313, 3233, 0, 1)}, sources, func);

//     // double le = (*mod)("YL",10*(gen-1)+gen-1);
//     // double G3=thdm_params->ld*(thdm_params->ld*thdm_params->lu+1.)*F6SP(thdm_params->xt,thdm_params->xH)+thdm_params->ld*thdm_params->lu*thdm_params->lu*F7SP(thdm_params->xt,thdm_params->xH)
// 	// +thdm_params->lu*thdm_params->lu*(thdm_params->ld*F8SP(thdm_params->xt,thdm_params->xH)+thdm_params->lu*F9SP(thdm_params->xt,thdm_params->xH)+thdm_params->lu*F10SP(thdm_params->xt,thdm_params->xH))+thdm_params->lu*F11SP(thdm_params->xt,thdm_params->xH)+thdm_params->lu*F12SP(thdm_params->xt,thdm_params->xH);
//     // double CPn_2HDM=thdm_params->xt*(-le*(thdm_params->ld*F1SP(thdm_params->xt,thdm_params->xH)+thdm_params->lu*F2SP(thdm_params->xt,thdm_params->xH))+le*thdm_params->lu*F3SP(thdm_params->xt,thdm_params->xH))+thdm_params->xt/2./thdm_params->xA*(le)*G3;

//     // double coeff_temp=CPc_2HDM(thdm_params->xH,thdm_params->xt,thdm_params->lu,thdm_params->ld,le,thdm_params->sw2)+CPn_2HDM;
//     // coeff_temp*=(wilson_p("WPARAM_SI_SM", 3)*thdm_params->mass_b_muW/sm("MASS",24)/sm("MASS",24))/thdm_params->sw2;
//     // // return this->double_to_complex_save("LO", coeff_temp);
// }

// void C_Blnu_P_THDM::LO_calculation() {

//     std::unordered_set<ParamId> sources {
//         {"WPARAM_SI_SM", 4},
//         {"WPARAM_SI_BSM", 7},
//         {"WPARAM_SI_BSM", 8},
//         {ParameterType::SM, "MASS", 15},
//         {ParameterType::BSM, "MASS", 37},
//         {ParameterType::BSM, "YL", {2, 2}},
//         {ParameterType::SM, "QCD", {5, 1}},
//     };

//     auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
//         double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
//         double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
//         double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();

//         double m_b = src.at({ParameterType::SM, "QCD", {5, 1}})->get_val();
//         double m_tau = src.at({ParameterType::SM, "MASS", 15})->get_val();
//         double l_tau = src.at({ParameterType::BSM, "YL", {2, 2}})->get_val();

//         dep_param->set_expected(-m_b * m_tau * ld * l_tau / std::pow(mH, 2));
//     };

//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(2051516, 3434, 0, 1)}, sources, func);

//     // double m_b = QCDHelper::mass_b_msbar();
//     // double m_tau = sm("MASS", 15);
//     // double l_tau = (*mod)("YL", 22);
//     // // return this->double_to_complex_save("LO", -m_b * m_tau * thdm_params->ld * l_tau / std::pow(thdm_params->m_H, 2));
// }

// // TODO : need to merge B_lnu group and B_CLNU group
// void C_S1_THDM::LO_calculation() {

//     std::unordered_set<ParamId> sources {
//         {"WPARAM_SI_SM", 4},
//         {"WPARAM_SI_BSM", 7},
//         {"WPARAM_SI_BSM", 8},
//         {ParameterType::SM, "MASS", 15},
//         {ParameterType::BSM, "MASS", 37},
//         {ParameterType::BSM, "YL", {2, 2}},
//         {ParameterType::SM, "QCD", {5, 1}},
//     };


//     auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
//         double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
//         double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
//         double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();

//         double m_b = (*Parameters::GetInstance())("QCD", LhaID(5, 1));
//         double m_tau = src.at({ParameterType::SM, "MASS", 15})->get_val();
//         double l_tau = src.at({ParameterType::BSM, "YL", 22})->get_val();

//         dep_param->set_expected(-m_b * m_tau * ld * l_tau / std::pow(mH, 2));
//     };

//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(4051516, 3231, 0, 1)}, sources, func);

//     // double m_b = QCDHelper::mass_b_msbar();
//     // double m_tau = sm("MASS", 15);
//     // double l_tau = (*mod)("YL", 22);
//     // // return this->double_to_complex_save("LO", -m_b * m_tau * thdm_params->ld * l_tau / std::pow(thdm_params->m_H, 2));
// }

// void C_S2_THDM::LO_calculation() {

//     std::unordered_set<ParamId> sources {
//         {"WPARAM_SI_SM", 4},
//         {"WPARAM_SI_BSM", 7},
//         {"WPARAM_SI_BSM", 8},
//         {ParameterType::SM, "MASS", 15},
//         {ParameterType::BSM, "MASS", 37},
//         {ParameterType::BSM, "YL", {2, 2}},
//         {ParameterType::SM, "MASS", 4},
//     };

//     auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
//         double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
//         double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
//         double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();

//         double m_c = src.at({ParameterType::SM, "MASS", 4})->get_val(); // TODO : mass c does not run ?
//         double m_tau = src.at({ParameterType::SM, "MASS", 15})->get_val();
//         double l_tau = src.at({ParameterType::BSM, "YL", {2, 2}})->get_val();

//         dep_param->set_expected(-m_c * m_tau * lu * l_tau / std::pow(mH, 2));
//     };

//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(4051516, 3131, 0, 1)}, sources, func);

//     // double m_c = sm("MASS", 4);
//     // double m_tau = sm("MASS", 15);
//     // double l_tau = (*mod)("YL", 22);
//     // // return this->double_to_complex_save("LO", -m_c * m_tau * thdm_params->lu * l_tau / std::pow(thdm_params->m_H, 2));
// }

// void WilsonCoefficient_THDM::init(QCDOrder order) {
//     if (!is_owned) {
//         thdm_parameters::init();
//     }

//     WilsonCoefficient::init(order);
// }

#include "MartyWilsonSuper.h"

// void MartyWilson::LO_calculation() {

//     LOG_INFO("In MartyWilson::LO_calculation of coefficient " + this->get_name());

//     std::unordered_set<ParamId> sources {
//         {"EW_SCALE", 1}
//     };

//     auto func = [this, &sources] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         LOG_INFO("Updating coeff");
//         double epsi = 1e-4;
//         double ew_scale = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
//         complex_t result = 0;

//         // for (size_t i = 0; i < df.getRowCount(); ++i) {
//         //     double Q_match = df.iat<double>(i, "Q_match");
//         //     if (fabs(Q_match - ew_scale) < epsi) {
//         //         std::cout << this->get_name() << " waw" << std::endl;
//         //         for (auto& _ : this->df.getColumnNames()) {
//         //             if (this->get_name()+"_real" == _) {
//         //                 if (isnan(df.iat<double>(i, this->get_name()+"_real")) && isnan(df.iat<double>(i, this->get_name()+"_img"))) {
//         //                     break;
//         //                 }
//         //                 std::cout << df.iat<double>(i, this->get_name()+"_real") << " BUTE" << std::endl;
//         //                 result = {df.iat<double>(i, this->get_name()+"_real"), df.iat<double>(i, this->get_name()+"_img")};
//         //                 dep_param->set_expected(result);
//         //                 return;
//         //                 // this->set_WilsonCoeffMatching("LO", {df.iat<double>(i, this->get_name()+"_real"), df.iat<double>(i, this->get_name()+"_img")});
//         //                 //return {df.iat<double>(i, this->get_name()+"_real"), df.iat<double>(i, this->get_name()+"_img")}; TODO
//         //             }
//         //         } 
//         //     }
//         // }

//         MartyInterface martyInterface;
//         martyInterface.calculate(this->get_name(), this->get_model(), ew_scale);
//         df = csv_reader.read_csv(this->csv_path);
//         df.setIndex(df.getColumn<double>("Q_match").to_string_vec());

//         for (size_t i = 0; i < df.getRowCount(); ++i) {
//             double Q_match = df.iat<double>(i, "Q_match");
//             if (fabs(Q_match - ew_scale) < epsi) {
//                 result = {df.iat<double>(i, this->get_name()+"_real"), df.iat<double>(i, this->get_name()+"_img")};
//                 break;
//                 // this->set_WilsonCoeffMatching("LO", {df.iat<double>(i, this->get_name()+"_real"), df.iat<double>(i, this->get_name()+"_img")});
//                 //return {df.iat<double>(i, this->get_name()+"_real"), df.iat<double>(i, this->get_name()+"_img")}; 
//             }
//         }

//         std::set<std::string> special = {"KIN", "WEIN", "Finite", "REGPROP"}; //TODO : do better with this, in SMParamSetter
//         for (auto &par : martyInterface.get_dependencies(this->get_name())) {
//             if (std::find(special.begin(), special.end(),par.block) != special.end()) {
//                 continue;
//             }
//             if (par.is_bsm) {
//                 sources.emplace(ParamId{ParameterType::BSM, par.block, par.code});
//             } else {
//                 sources.emplace(ParamId{ParameterType::SM, par.block, par.code});
//             }
//         }

//         dep_param->set_expected(result);
//     };
//     ParamId pid {ParameterType::WILSON, "EW_SCALE", 1};
//     std::unordered_map<ParamId, std::shared_ptr<Parameter>> dummy {{pid, std::make_shared<Parameter>(pid, 1, 0, 0)}};
//     func(dummy, std::make_shared<DependentParameter>(dummy, func));

//     //TODO  :configuration here to SM, need to put at FULL but cannot work with WilsonGroup at the moment (only sm)
//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, WCoefMapper::flha_full(WCoefMapper::enum_elt(this->coeffName), QCDOrder::LO, this->type)}, sources, func);

// }
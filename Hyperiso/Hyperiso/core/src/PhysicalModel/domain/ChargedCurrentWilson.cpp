// #include "ChargedCurrentWilson.h"

// void C_Blnu_A::LO_calculation() {
//     std::unordered_set<ParamId> sources {};

//     auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         dep_param->set_expected(1.);
//     };

//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, WCoefMapper::flha_full(WCoef::CBlnu_A, QCDOrder::LO, this->type)}, sources, func);
// }

// void C_V1::LO_calculation() {
//     std::unordered_set<ParamId> sources {};

//     auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
//         dep_param->set_expected(1.);
//     };

//     WilsonParamComposer().compose_parameter(ParamId{this->storage_block, WCoefMapper::flha_full(WCoef::C_V1, QCDOrder::LO, this->type)}, sources, func);
// }

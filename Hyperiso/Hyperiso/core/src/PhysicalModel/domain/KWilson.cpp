#include "KWilson.h"

CK9::CK9() : WilsonCoefficient("CK9", GroupMapper::str(WGroup::K) + "_MATCH") {
    // LO
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)}  // x_t
        },
        compute_LO,
        LhaID(0, 1, 0, 0)
    };

}

double CK9::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
    return -0.5 * F0t(xt) - 1. / 3.;
}

CK10::CK10() : WilsonCoefficient("CK10", GroupMapper::str(WGroup::K) + "_MATCH") {
    // LO
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)}  // x_t
        },
        compute_LO,
        LhaID(0, 2, 0, 0)
    };

}

double CK10::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
    return -0.5 * F0t(xt) - 1. / 3.;
}

CKQ1::CKQ1() : WilsonCoefficient("CKQ1", GroupMapper::str(WGroup::K) + "_MATCH") {
    // LO
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)}  // x_t
        },
        compute_LO,
        LhaID(0, 3, 0, 0)
    };

}

double CKQ1::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
    return -0.5 * F0t(xt) - 1. / 3.;
}

CKQ2::CKQ2() : WilsonCoefficient("CKQ2", GroupMapper::str(WGroup::K) + "_MATCH") {
    // LO
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)}  // x_t
        },
        compute_LO,
        LhaID(0, 4, 0, 0)
    };

}

double CKQ2::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
    return -0.5 * F0t(xt) - 1. / 3.;
}

CPK9::CPK9() : WilsonCoefficient("CPK9", GroupMapper::str(WGroup::K) + "_MATCH") {
    // LO
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)}  // x_t
        },
        compute_LO,
        LhaID(0, 5, 0, 0)
    };

}

double CPK9::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
    return -0.5 * F0t(xt) - 1. / 3.;
}

CPK10::CPK10() : WilsonCoefficient("CPK10", GroupMapper::str(WGroup::K) + "_MATCH") {
    // LO
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)}  // x_t
        },
        compute_LO,
        LhaID(0, 6, 0, 0)
    };

}

double CPK10::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
    return -0.5 * F0t(xt) - 1. / 3.;
}

CPKQ1::CPKQ1() : WilsonCoefficient("CPKQ1", GroupMapper::str(WGroup::K) + "_MATCH") {
    // LO
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)}  // x_t
        },
        compute_LO,
        LhaID(0, 7, 0, 0)
    };

}

double CPKQ1::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
    return -0.5 * F0t(xt) - 1. / 3.;
}

CPKQ2::CPKQ2() : WilsonCoefficient("CPKQ2", GroupMapper::str(WGroup::K) + "_MATCH") {
    // LO
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)}  // x_t
        },
        compute_LO,
        LhaID(0, 8, 0, 0)
    };

}

double CPKQ2::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
    return -0.5 * F0t(xt) - 1. / 3.;
}
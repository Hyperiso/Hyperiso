#include "KWilson.h"

CK9::CK9() : WilsonCoefficient("CK9", GroupMapper::str(WGroup::K, ScaleType::MATCHING)) {
    // LO
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)}  // x_t
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };

}

double CK9::compute_LO(const ParamSrc& src) {
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    return -0.5 * F0t(xt) - 1. / 3.;
}

CK10::CK10() : WilsonCoefficient("CK10", GroupMapper::str(WGroup::K, ScaleType::MATCHING)) {
    // LO
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)},  // x_t
            {"WPARAM_SI_SM", 4}                // sw2
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };

    matching_info[QCDOrder::NLO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)},  // x_t
            {"WPARAM_SI_SM", 4},               // sw2
            {ParameterType::SM, "MASS", 24},   // mW
            {ParameterType::SM, "QCD", 6},     //mt 
        },
        compute_NLO,
        get_lhaid_from_name(QCDOrder::NLO)
    };

}

double CK10::compute_LO(const ParamSrc& src) {
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    double sw2  = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);
    return -Y0(xt)/sw2;
}

//TODO : careful with NLO, alpha/4pi is not the same than before (not same scale)
double CK10::compute_NLO(const ParamSrc& src) {
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    double sw2  = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);
    double mW    = src.get_val(ParameterType::SM, "MASS", 24);
    double mtmt = src.get_val(ParameterType::SM, "QCD", 6);
    return -Y1(xt,mtmt,mW)/sw2;
}

CKQ1::CKQ1() : WilsonCoefficient("CKQ1", GroupMapper::str(WGroup::K, ScaleType::MATCHING)) {
    // LO
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)}  // x_t
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };

}

double CKQ1::compute_LO(const ParamSrc& src) {
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    return -0.5 * F0t(xt) - 1. / 3.;
}

CKQ2::CKQ2() : WilsonCoefficient("CKQ2", GroupMapper::str(WGroup::K, ScaleType::MATCHING)) {
    // LO
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)}  // x_t
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };

}

double CKQ2::compute_LO(const ParamSrc& src) {
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    return -0.5 * F0t(xt) - 1. / 3.;
}

CPK9::CPK9() : WilsonCoefficient("CPK9", GroupMapper::str(WGroup::K, ScaleType::MATCHING)) {
    // LO
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)}  // x_t
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };

}

double CPK9::compute_LO(const ParamSrc& src) {
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    return -0.5 * F0t(xt) - 1. / 3.;
}

CPK10::CPK10() : WilsonCoefficient("CPK10", GroupMapper::str(WGroup::K, ScaleType::MATCHING)) {
    // LO
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)}  // x_t
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };

}

double CPK10::compute_LO(const ParamSrc& src) {
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    return -0.5 * F0t(xt) - 1. / 3.;
}

CPKQ1::CPKQ1() : WilsonCoefficient("CPKQ1", GroupMapper::str(WGroup::K, ScaleType::MATCHING)) {
    // LO
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)}  // x_t
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };

}

double CPKQ1::compute_LO(const ParamSrc& src) {
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    return -0.5 * F0t(xt) - 1. / 3.;
}

CPKQ2::CPKQ2() : WilsonCoefficient("CPKQ2", GroupMapper::str(WGroup::K, ScaleType::MATCHING)) {
    // LO
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)}  // x_t
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };

}

double CPKQ2::compute_LO(const ParamSrc& src) {
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    return -0.5 * F0t(xt) - 1. / 3.;
}

CK_L::CK_L() : WilsonCoefficient("CK_L", GroupMapper::str(WGroup::K, ScaleType::MATCHING)) {
    // LO
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)},  // x_t
            {"WPARAM_SI_SM", 4},               // sw2
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };

    matching_info[QCDOrder::NLO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)},  // x_t
            {"WPARAM_SI_SM", 4},               // sw2
            {ParameterType::SM, "MASS", 24},   // mW
            {ParameterType::SM, "QCD", 6},     //mt 
        },
        compute_NLO,
        get_lhaid_from_name(QCDOrder::NLO)
    };

}

double CK_L::compute_LO(const ParamSrc& src) {
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    double sw2  = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);
    return X0(xt)/std::pow(sw2,2);
}

double CK_L::compute_NLO(const ParamSrc& src) {
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    double sw2  = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);
    double mW    = src.get_val(ParameterType::SM, "MASS", 24);
    double mtmt = src.get_val(ParameterType::SM, "QCD", 6);
    return X1(xt, mtmt, mW)/std::pow(sw2,2);
}
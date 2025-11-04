#include "BPrimeWilson.h"


// ---------- C'7 ----------

CP7::CP7() : WilsonCoefficient("CP7", GroupMapper::str(WGroup::BPrime, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", {2, 1}},     // x_t
            {"WPARAM_MATCH_SM", {5, 1}},     // m_b(muW)
            {ParameterType::SM, "MASS", 3}   // m_s
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };

    matching_info[QCDOrder::NLO] = MatchingInfo(get_lhaid_from_name(QCDOrder::NLO));

    matching_info[QCDOrder::NNLO] = MatchingInfo(get_lhaid_from_name(QCDOrder::NNLO));
}

double CP7::compute_LO(const ParamSrc& src) {
    double xt     = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});;
    double mb     = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1});;
    double ms     = src.get_val(ParameterType::SM, "MASS", 3);;

    printf("ms in SM (LO) : %.8lf\n", ms);
    printf("mb in SM (LO) : %.8lf\n", mb);
    printf("xt : %.8lf\n", xt);

    printf("CP7 in SM (LO) : %.8lf\n", ms / mb * (-0.5 * A0t(xt) - 23. / 36.));
    return ms / mb * (-0.5 * A0t(xt) - 23. / 36.);
}

// ---------- C'8 ----------

CP8::CP8() : WilsonCoefficient("CP8", GroupMapper::str(WGroup::BPrime, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", {2, 1}},     // x_t
            {"WPARAM_MATCH_SM", {5, 1}},     // m_b(muW)
            {ParameterType::SM, "MASS", 3}   // m_s
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };

    matching_info[QCDOrder::NLO] = MatchingInfo(get_lhaid_from_name(QCDOrder::NLO));

    matching_info[QCDOrder::NNLO] = MatchingInfo(get_lhaid_from_name(QCDOrder::NNLO));
}

double CP8::compute_LO(const ParamSrc& src) {
    double xt     = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});;
    double mb     = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1});;
    double ms     = src.get_val(ParameterType::SM, "MASS", 3);;
    printf("CP8 in SM (LO) : %.8lf\n", ms / mb * (-0.5 * F0t(xt) - 1. / 3.));
    return ms / mb * (-0.5 * F0t(xt) - 1. / 3.);
}
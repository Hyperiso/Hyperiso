#include "BScalarWilson.h"

// ---------- CQ1 ----------

CQ1::CQ1() : WilsonCoefficient("CQ1", GroupMapper::str(WGroup::BScalar, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", {2, 1}},     // x_t
            {"WPARAM_MATCH_SM", {5, 1}},     // m_b(muW) ?! (voir TODO)
            {"WPARAM_SI_SM", 1},             // xh
            {"WPARAM_SI_SM", 3},             // ml
            {"WPARAM_SI_SM", 4},             // sw2
            {ParameterType::SM, "MASS", 24}  // m_W
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };

    matching_info[QCDOrder::NLO] = MatchingInfo(get_lhaid_from_name(QCDOrder::NLO));

    matching_info[QCDOrder::NNLO] = MatchingInfo(get_lhaid_from_name(QCDOrder::NNLO));
}

double CQ1::compute_LO(const ParamSrc& src) {
    double xt     = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});;
    double mb_muW = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1});; // TODO : Ask Nazila (check SI first) : Why {5,2} and not {5,1} ?
    double xh     = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 1);;
    double ml     = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 3);;
    double sw2    = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);;
    double mW     = src.get_val(ParameterType::SM, "MASS", 24);;

    // printf("\n----- In CQ1 -----\n");
    // printf("xt = %.4e\n", xt);
    // printf("mb_muW = %.4e\n", mb_muW);
    // printf("xh = %.4e\n", xh);
    // printf("ml = %.4e\n", ml);
    // printf("sw2 = %.4e\n", sw2);
    // printf("mW = %.4e\n", mW);

    double CSc_SM = -xt * (xt - 2.) / 12. / pow(xt - 1., 2)
                  + (xt - 2.) * (3. * xt - 1.) / 24. / pow(xt - 1., 3.) * log(xt);

    double CSn_SMonly = -3. * xt / 8. / xh + xt * F0SP(xt);

    double coeff_temp = (CSc_SM + CSn_SMonly) * (ml * mb_muW / (mW * mW)) / sw2;

    // printf("CSc = %.4e\n", CSc_SM);
    // printf("CSn = %.4e\n", CSn_SMonly);
    // printf("CQ1 = %.4e\n", coeff_temp);

    return coeff_temp;
}

// ---------- CQ2 ----------

CQ2::CQ2() : WilsonCoefficient("CQ2", GroupMapper::str(WGroup::BScalar, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", {2, 1}},  // xt
            {"WPARAM_MATCH_SM", {2, 2}},  // xt^2
            {"WPARAM_MATCH_SM", {2, 3}},  // xt^3
            {"WPARAM_MATCH_SM", {2, 4}},  // xt^4
            {"WPARAM_MATCH_SM", {5, 1}},  // m_b(muW)
            {"WPARAM_SI_SM", 3},          // ml
            {"WPARAM_SI_SM", 4},          // sw2
            {ParameterType::SM, "MASS", 24} // m_W
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };

    matching_info[QCDOrder::NLO] = MatchingInfo(get_lhaid_from_name(QCDOrder::NLO));

    matching_info[QCDOrder::NNLO] = MatchingInfo(get_lhaid_from_name(QCDOrder::NNLO));
}

double CQ2::compute_LO(const ParamSrc& src) {
    double xt   = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    double xt2  = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 2});
    double xt3  = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 3});
    double xt4  = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 4});
    double mb   = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1});
    double ml   = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 3);
    double sw2  = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);
    double mW   = src.get_val(ParameterType::SM, "MASS", 24);

    // printf("\n----- In CQ2 -----\n");
    // printf("xt = %.4e\n", xt);
    // printf("mb = %.4e\n", mb);
    // printf("ml = %.4e\n", ml);
    // printf("sw2 = %.4e\n", sw2);
    // printf("mW = %.4e\n", mW);

    double CPc_SM =
        1. / 24. * (
            xt * (36. * xt3 - 203. * xt2 + 352. * xt - 209.) / 6. / pow(xt - 1., 3.)
            + (17. * xt4 - 34. * xt3 + 4. * xt2 + 23. * xt - 6.) / pow(xt - 1., 4.) * log(xt)
        )
        - sw2 / 36. * (
            xt * (18. * xt3 - 139. * xt2 + 274. * xt - 129.) / 2. / pow(xt - 1., 3.)
            + (24. * xt4 - 33. * xt3 - 45. * xt2 + 50. * xt - 8.) / pow(xt - 1., 4.) * log(xt)
        );

    double CPn_SMonly = 0.;

    // printf("CPc = %.4e\n", CPc_SM);
    // printf("CPn = %.4e\n", CPn_SMonly);
    // printf("CQ2 = %.4e\n", (CPc_SM + CPn_SMonly) * (ml * mb / (mW * mW)) / sw2);

    return (CPc_SM + CPn_SMonly) * (ml * mb / (mW * mW)) / sw2;
}
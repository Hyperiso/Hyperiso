#include "BWilsonTHDM.h"

C3_THDM::C3_THDM() : WilsonCoefficient("C3_THDM", GroupMapper::str(WGroup::B) + "_MATCH") {
    matching_info[QCDOrder::NNLO] = {
        {
            {"WPARAM_SI_BSM", 7},            // lu
            {"WPARAM_MATCH_BSM", 1},         // yt
            {"EW_SCALE", 1},                 // Q_match
            {ParameterType::BSM, "MASS", 37} // m_H
        },
        compute_NNLO,
        LhaID(3050707, 4133, 2, 1)           // LHA ID spécifique au THDM
    };
}

double C3_THDM::compute_NNLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double lu     = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
    double yt     = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
    double Q      = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
    double mH     = src.at({ParameterType::BSM, "MASS", 37})->get_val();

    return G3H(yt, lu) + Delta3H(yt, lu) * log(pow(Q / mH, 2.0));
}

C4_THDM::C4_THDM() : WilsonCoefficient("C4_THDM", GroupMapper::str(WGroup::B) + "_MATCH") {
    // NLO
    matching_info[QCDOrder::NLO] = {
        {
            {"WPARAM_SI_BSM", 7},           // lu
            {"WPARAM_MATCH_BSM", 1}        // yt
        },
        compute_NLO,
        LhaID(3050707, 6153, 1, 1)
    };

    // NNLO
    matching_info[QCDOrder::NNLO] = {
        {
            {"WPARAM_SI_BSM", 7},            // lu
            {"WPARAM_MATCH_BSM", 1},         // yt
            {"EW_SCALE", 1},                 // Q_match
            {ParameterType::BSM, "MASS", 37} // m_H
        },
        compute_NNLO,
        LhaID(3050707, 6153, 2, 1)
    };
}

double C4_THDM::compute_NLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
    double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();

    return EH(yt, lu);
}

double C4_THDM::compute_NNLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double lu     = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
    double yt     = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
    double Q      = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
    double mH     = src.at({ParameterType::BSM, "MASS", 37})->get_val();

    return G4H(yt, lu) + Delta4H(yt, lu) * log(pow(Q / mH, 2.0));
}

C5_THDM::C5_THDM() : WilsonCoefficient("C5_THDM", GroupMapper::str(WGroup::B) + "_MATCH") {
    matching_info[QCDOrder::NNLO] = {
        {
            {"WPARAM_SI_BSM", 7},            // lu
            {"WPARAM_MATCH_BSM", 1},         // yt
            {"EW_SCALE", 1},                 // Q_match
            {ParameterType::BSM, "MASS", 37} // m_H
        },
        compute_NNLO,
        LhaID(3050707, 4536, 2, 1)
    };
}

double C5_THDM::compute_NNLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double lu     = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
    double yt     = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
    double Q      = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
    double mH     = src.at({ParameterType::BSM, "MASS", 37})->get_val();

    double C4H_1 = EH(yt, lu);
    double C3H_2 = G3H(yt, lu) + Delta3H(yt, lu) * log(pow(Q / mH, 2.0));

    return -C3H_2 / 10. + 2. / 15. * C4H_1;
}

C6_THDM::C6_THDM() : WilsonCoefficient("C6_THDM", GroupMapper::str(WGroup::B) + "_MATCH") {
    matching_info[QCDOrder::NNLO] = {
        {
            {"WPARAM_SI_BSM", 7},            // lu
            {"WPARAM_MATCH_BSM", 1},         // yt
            {"EW_SCALE", 1},                 // Q_match
            {ParameterType::BSM, "MASS", 37} // m_H
        },
        compute_NNLO,
        LhaID(3050707, 6556, 2, 1)
    };
}

double C6_THDM::compute_NNLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double lu     = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
    double yt     = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
    double Q      = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
    double mH     = src.at({ParameterType::BSM, "MASS", 37})->get_val();

    double C4H_1 = EH(yt, lu);
    double C3H_2 = G3H(yt, lu) + Delta3H(yt, lu) * log(pow(Q / mH, 2.0));

    return -3. / 16. * C3H_2 + 1. / 4. * C4H_1;
}

C7_THDM::C7_THDM() : WilsonCoefficient("C7_THDM", GroupMapper::str(WGroup::B) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_SI_BSM", 7},           // lu
            {"WPARAM_SI_BSM", 8},           // ld
            {"WPARAM_MATCH_BSM", 1}         // yt
        },
        compute_LO,
        LhaID(305, 4422, 0, 1)
    };

    matching_info[QCDOrder::NLO] = {
        {
            {"WPARAM_SI_BSM", 7},           // lu
            {"WPARAM_SI_BSM", 8},           // ld
            {"WPARAM_MATCH_BSM", 1},        // yt
            {"EW_SCALE", 1},                // Q_match
            {ParameterType::BSM, "MASS", 37} // m_H
        },
        compute_NLO,
        LhaID(305, 4422, 1, 1)
    };

    matching_info[QCDOrder::NNLO] = {
        {
            {"WPARAM_SI_BSM", 7},           // lu
            {"WPARAM_SI_BSM", 8},           // ld
            {"WPARAM_MATCH_BSM", 1},        // yt
            {"WPARAM_MATCH_SM", 6},         // mass_top_muW
            {"EW_SCALE", 1}                 // Q_match
        },
        compute_NNLO,
        LhaID(305, 4422, 2, 1)
    };
}

double C7_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
    double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
    double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();

    return 1. / 3. * lu * lu * F7_1(yt) - lu * ld * F7_2(yt);
}

double C7_THDM::compute_NLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
    double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
    double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
    double Q  = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
    double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();

    return G7H(yt, lu, ld)
         + Delta7H(yt, lu, ld) * log(pow(Q / mH, 2.0))
         - 4. / 9. * EH(yt, lu);
}

double C7_THDM::compute_NNLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double lu    = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
    double ld    = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
    double yt    = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
    double mtop  = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
    double Q     = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
    std::cout << "mtop : " << mtop << std::endl;
    return C7H2(yt, lu, ld, log(pow(Q / mtop, 2.0)));
}

C8_THDM::C8_THDM() : WilsonCoefficient("C8_THDM", GroupMapper::str(WGroup::B) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_SI_BSM", 7},           // lu
            {"WPARAM_SI_BSM", 8},           // ld
            {"WPARAM_MATCH_BSM", 1}         // yt
        },
        compute_LO,
        LhaID(305, 6421, 0, 1)
    };

    matching_info[QCDOrder::NLO] = {
        {
            {"WPARAM_SI_BSM", 7},           // lu
            {"WPARAM_SI_BSM", 8},           // ld
            {"WPARAM_MATCH_BSM", 1},        // yt
            {"EW_SCALE", 1},                // Q_match
            {ParameterType::BSM, "MASS", 37} // m_H
        },
        compute_NLO,
        LhaID(305, 6421, 1, 1)
    };

    matching_info[QCDOrder::NNLO] = {
        {
            {"WPARAM_SI_BSM", 7},           // lu
            {"WPARAM_SI_BSM", 8},           // ld
            {"WPARAM_MATCH_BSM", 1},        // yt
            {"WPARAM_MATCH_SM", 6},         // mtop(muW)
            {"EW_SCALE", 1}                 // Q_match
        },
        compute_NNLO,
        LhaID(305, 6421, 2, 1)
    };
}

double C8_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
    double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
    double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();

    return 1. / 3. * lu * ld * F8_1(yt) - lu * ld * F8_2(yt);
}

double C8_THDM::compute_NLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
    double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
    double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
    double Q  = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
    double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();

    return G8H(yt, lu, ld)
         + Delta8H(yt, lu, ld) * log(pow(Q / mH, 2.0))
         - 1. / 6. * EH(yt, lu);
}

double C8_THDM::compute_NNLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double lu    = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
    double ld    = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
    double yt    = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
    double mtop  = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
    double Q     = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();

    return C8H2(yt, lu, ld, log(pow(Q / mtop, 2.0)));
}

C9_THDM::C9_THDM() : WilsonCoefficient("C9_THDM", GroupMapper::str(WGroup::B) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_SI_SM", 4},               // sw2
            {"WPARAM_SI_BSM", 7},              // lu
            {"WPARAM_MATCH_BSM", 1},           // yt
            {"WPARAM_MATCH_SM", {2, 1}}        // xt
        },
        compute_LO,
        LhaID(3051313, 4133, 0, 1)
    };

    matching_info[QCDOrder::NLO] = {
        {
            {"WPARAM_SI_SM", 4},               // sw2
            {"WPARAM_SI_BSM", 7},              // lu
            {"WPARAM_MATCH_BSM", 1},           // yt
            {"WPARAM_MATCH_SM", {2, 1}},       // xt
            {"EW_SCALE", 1},                   // Q_match
            {ParameterType::BSM, "MASS", 37}   // m_H
        },
        compute_NLO,
        LhaID(3051313, 4133, 1, 1)
    };
}


double C9_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
    double lu  = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
    double yt  = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
    double xt  = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();

    std::cout << "sw2 : " << sw2 << std::endl;
    std::cout << "xt : " << xt << std::endl;
    std::cout << "yt : " << yt << std::endl;
    return (1. - 4. * sw2) / sw2 * C9llH0(xt, yt, lu) - D9H0(yt, lu);
}

double C9_THDM::compute_NLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
    double lu  = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
    double yt  = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
    double xt  = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
    double Q   = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
    double mH  = src.at({ParameterType::BSM, "MASS", 37})->get_val();

    double logQH2 = log(pow(Q / mH, 2.0));

    return (1. - 4. * sw2) / sw2 * C9llH1(xt, yt, lu, logQH2)
         - D9H1(yt, lu, logQH2);
}

C10_THDM::C10_THDM() : WilsonCoefficient("C10_THDM", GroupMapper::str(WGroup::B) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_SI_SM", 4},               // sw2
            {"WPARAM_SI_BSM", 7},              // lu
            {"WPARAM_MATCH_BSM", 1},           // yt
            {"WPARAM_MATCH_SM", {2, 1}}        // xt
        },
        compute_LO,
        LhaID(3051313, 4137, 0, 1)
    };

    matching_info[QCDOrder::NLO] = {
        {
            {"WPARAM_SI_SM", 4},               // sw2
            {"WPARAM_SI_BSM", 7},              // lu
            {"WPARAM_MATCH_BSM", 1},           // yt
            {"WPARAM_MATCH_SM", {2, 1}},       // xt
            {"EW_SCALE", 1},                   // Q_match
            {ParameterType::BSM, "MASS", 37}   // m_H
        },
        compute_NLO,
        LhaID(3051313, 4137, 1, 1)
    };
}

double C10_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
    double lu  = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
    double yt  = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
    double xt  = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();

    return -C9llH0(xt, yt, lu) / sw2;
}

double C10_THDM::compute_NLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
    double lu  = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
    double yt  = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
    double xt  = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
    double Q   = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
    double mH  = src.at({ParameterType::BSM, "MASS", 37})->get_val();

    return -C9llH1(xt, yt, lu, log(pow(Q / mH, 2.0))) / sw2;
}

CQ1_THDM::CQ1_THDM() : WilsonCoefficient("CQ1_THDM", GroupMapper::str(WGroup::B) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_SI_SM", 3},              // ml
            {"WPARAM_SI_SM", 4},              // sw2
            {"WPARAM_SI_BSM", 1},             // xh
            {"WPARAM_SI_BSM", 2},             // xH
            {"WPARAM_SI_BSM", 3},             // xH0
            {"WPARAM_SI_BSM", 6},             // beta
            {"WPARAM_SI_BSM", 7},             // lu
            {"WPARAM_SI_BSM", 8},             // ld
            {"WPARAM_SI_BSM", 9},             // alpha
            {"WPARAM_SI_BSM", 10},            // le
            {"WPARAM_MATCH_SM", {2, 1}},      // xt
            {"WPARAM_MATCH_SM", {5, 1}},      // mass_b_muW
            {ParameterType::SM, "MASS", 24}   // m_W
        },
        compute_LO,
        LhaID(3051313, 3230, 0, 1)
    };
}


double CQ1_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double xh     = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 1})->get_val();
    double ml     = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 3})->get_val();
    double sw2    = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
    double xH     = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 2})->get_val();
    double xH0    = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 3})->get_val();
    double beta   = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val();
    double lu     = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
    double ld     = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
    double alpha  = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 9})->get_val();
    double le     = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 10})->get_val();
    double xt     = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
    double mb_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
    double mW     = src.at({ParameterType::SM, "MASS", 24})->get_val();

    double G1 = -3. / 4. + ld * lu * F4SP(xt, xH) + lu * lu * F5SP(xt, xH);
    double G2 = ld * (ld * lu + 1.) * F6SP(xt, xH)
              - ld * lu * lu * F7SP(xt, xH)
              + lu * lu * (ld * F8SP(xt, xH) + lu * F9SP(xt, xH) - lu * F10SP(xt, xH))
              + lu * F11SP(xt, xH) - lu * F12SP(xt, xH);

    double CSn_2HDM =
        xt * (F0SP(xt) + le * (ld * F1SP(xt, xH) + lu * F2SP(xt, xH)) + le * lu * F3SP(xt, xH)) +
        xt / (2. * xh) * (sin(alpha - beta) + cos(alpha - beta) * le) *
        (sin(alpha - beta) * G1 + cos(alpha - beta) * G2) +
        xt / (2. * xH0) * (cos(alpha - beta) - sin(alpha - beta) * le) *
        (cos(alpha - beta) * G1 - sin(alpha - beta) * G2);

        printf("first_part : %.9lf\n", xt * (F0SP(xt) + le * (ld * F1SP(xt, xH) + lu * F2SP(xt, xH)) + le * lu * F3SP(xt, xH)));
        printf("second_part : %.9lf\n", xt / (2. * xh) * (sin(alpha - beta) + cos(alpha - beta) * le) *
        (sin(alpha - beta) * G1 + cos(alpha - beta) * G2));

    LOG_INFO("F0SP =", F0SP(xt));
    LOG_INFO("F1SP =", F1SP(xt, xH));
    LOG_INFO("F2SP =", F2SP(xt, xH));
    LOG_INFO("F3SP =", F3SP(xt, xH));
    LOG_INFO("alpha =", alpha);
    LOG_INFO("beta =", beta);
    LOG_INFO("xH0 =", xH0);
    LOG_INFO("xt =", xt);
    LOG_INFO("xh =", xh);
    LOG_INFO("xH =", xH);
    LOG_INFO("G1 =", G1);
    LOG_INFO("G2 =", G2);
    LOG_INFO("lu =", lu);
    LOG_INFO("ld =", ld);
    LOG_INFO("lambda_e =", le);
    LOG_INFO("sin(alpha-beta) =", sin(alpha - beta));
    LOG_INFO("cos(alpha-beta) =", cos(alpha - beta));
    LOG_INFO("beta =", beta);
    LOG_INFO("CSn_2HDM =", CSn_2HDM);
    LOG_INFO("CSc_2HDM =", CSc_2HDM(xH, xt, lu, ld, le));

    double coeff_temp = CSc_2HDM(xH, xt, lu, ld, le) + CSn_2HDM;
    coeff_temp *= (ml * mb_muW / (mW * mW)) / sw2;

    printf("sw2 = %lf\n", sw2);
    printf("ml = %lf\n", ml);
	printf("mass_b_muW = %.8lf\n", mb_muW);
	printf("param->mass_W = %lf\n", mW);
    std::cout << "CQ1H_0 : " << coeff_temp << std::endl;
    return coeff_temp;
}

CQ2_THDM::CQ2_THDM() : WilsonCoefficient("CQ2_THDM", GroupMapper::str(WGroup::B) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_SI_SM", 3},               // ml
            {"WPARAM_SI_SM", 4},               // sw2
            {"WPARAM_SI_BSM", 2},              // xH
            {"WPARAM_SI_BSM", 4},              // xA
            {"WPARAM_SI_BSM", 7},              // lu
            {"WPARAM_SI_BSM", 8},              // ld
            {"WPARAM_SI_BSM", 10},             // le
            {"WPARAM_MATCH_SM", {2, 1}},       // xt
            {"WPARAM_MATCH_SM", {5, 1}},       // mass_b_muW
            {ParameterType::SM, "MASS", 24}    // m_W
        },
        compute_LO,
        LhaID(3051313, 3233, 0, 1)
    };
}

double CQ2_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double ml       = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 3})->get_val();
    double sw2      = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
    double xH       = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 2})->get_val();
    double xA       = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 4})->get_val();
    double lu       = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
    double ld       = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
    double le       = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 10})->get_val();
    double xt       = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
    double mb_muW   = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
    double mW       = src.at({ParameterType::SM, "MASS", 24})->get_val();

    double G3 =
        ld * (ld * lu + 1.) * F6SP(xt, xH)
        + ld * lu * lu * F7SP(xt, xH)
        + lu * lu * (ld * F8SP(xt, xH) + lu * F9SP(xt, xH) + lu * F10SP(xt, xH))
        + lu * F11SP(xt, xH) + lu * F12SP(xt, xH);

    double CPn_2HDM =
        xt * (-le * (ld * F1SP(xt, xH) + lu * F2SP(xt, xH)) + le * lu * F3SP(xt, xH))
        + xt / (2. * xA) * le * G3;

    double coeff_temp = CPc_2HDM(xH, xt, lu, ld, le, sw2) + CPn_2HDM;
    coeff_temp *= (ml * mb_muW / (mW * mW)) / sw2;

    LOG_INFO("G3 =", G3);
    LOG_INFO("Factor =", (ml * mb_muW / (mW * mW)) / sw2);
    LOG_INFO("m_b_muW =", mb_muW);
    LOG_INFO("CPn_2HDM =", CPn_2HDM);
    LOG_INFO("CPc_2HDM =", CPc_2HDM(xH, xt, lu, ld, le, sw2));
    LOG_INFO("CQ2(mu_W) =", coeff_temp);

    return coeff_temp;
}


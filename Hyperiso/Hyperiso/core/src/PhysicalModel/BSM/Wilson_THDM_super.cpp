#include "Wilson_THDM_super.h"

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
            {"WPARAM_SI_BSM", 2},             // xH
            {"WPARAM_SI_BSM", 3},             // xH0
            {"WPARAM_SI_BSM", 6},             // beta
            {"WPARAM_SI_BSM", 7},             // lu
            {"WPARAM_SI_BSM", 8},             // ld
            {"WPARAM_SI_BSM", 9},             // alpha
            {"WPARAM_SI_BSM", 10},            // le
            {"WPARAM_MATCH_SM", {2, 1}},      // xt
            {"WPARAM_MATCH_SM", 6},           // mass_top_muW
            {ParameterType::SM, "MASS", 24}   // m_W
        },
        compute_LO,
        LhaID(3051313, 3230, 0, 1)
    };
}


double CQ1_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
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
    double mt_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
    double mW     = src.at({ParameterType::SM, "MASS", 24})->get_val();

    double G1 = -3. / 4. + ld * lu * F4SP(xt, xH) + lu * lu * F5SP(xt, xH);
    double G2 = ld * (ld * lu + 1.) * F6SP(xt, xH)
              - ld * lu * lu * F7SP(xt, xH)
              + lu * lu * (ld * F8SP(xt, xH) + lu * F9SP(xt, xH) - lu * F10SP(xt, xH))
              + lu * F11SP(xt, xH) - lu * F12SP(xt, xH);

    double CSn_2HDM =
        xt * (F0SP(xt) + le * (ld * F1SP(xt, xH) + lu * F2SP(xt, xH)) + le * lu * F3SP(xt, xH)) +
        xt / (2. * xH) * (sin(alpha - beta) + cos(alpha - beta) * le) *
        (sin(alpha - beta) * G1 + cos(alpha - beta) * G2) +
        xt / (2. * xH0) * (cos(alpha - beta) - sin(alpha - beta) * le) *
        (cos(alpha - beta) * G1 - sin(alpha - beta) * G2);

    double coeff_temp = CSc_2HDM(xH, xt, lu, ld, le) + CSn_2HDM;
    coeff_temp *= (ml * mt_muW / (mW * mW)) / sw2;

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

    return coeff_temp;
}


void C_Blnu_P_THDM::LO_calculation() {

    std::unordered_set<ParamId> sources {
        {"WPARAM_SI_SM", 4},
        {"WPARAM_SI_BSM", 7},
        {"WPARAM_SI_BSM", 8},
        {ParameterType::SM, "MASS", 15},
        {ParameterType::BSM, "MASS", 37},
        {ParameterType::BSM, "YL", {2, 2}},
        {ParameterType::SM, "QCD", {5, 1}},
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
        double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
        double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
        double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();

        double m_b = src.at({ParameterType::SM, "QCD", {5, 1}})->get_val();
        double m_tau = src.at({ParameterType::SM, "MASS", 15})->get_val();
        double l_tau = src.at({ParameterType::BSM, "YL", {2, 2}})->get_val();

        dep_param->set_expected(-m_b * m_tau * ld * l_tau / std::pow(mH, 2));
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(2051516, 3434, 0, 1)}, sources, func);

    // double m_b = QCDHelper::mass_b_msbar();
    // double m_tau = sm("MASS", 15);
    // double l_tau = (*mod)("YL", 22);
    // // return this->double_to_complex_save("LO", -m_b * m_tau * thdm_params->ld * l_tau / std::pow(thdm_params->m_H, 2));
}

// TODO : need to merge B_lnu group and B_CLNU group
void C_S1_THDM::LO_calculation() {

    std::unordered_set<ParamId> sources {
        {"WPARAM_SI_SM", 4},
        {"WPARAM_SI_BSM", 7},
        {"WPARAM_SI_BSM", 8},
        {ParameterType::SM, "MASS", 15},
        {ParameterType::BSM, "MASS", 37},
        {ParameterType::BSM, "YL", {2, 2}},
        {ParameterType::SM, "QCD", {5, 1}},
    };


    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
        double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
        double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
        double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();

        double m_b = (*Parameters::GetInstance())("QCD", LhaID(5, 1));
        double m_tau = src.at({ParameterType::SM, "MASS", 15})->get_val();
        double l_tau = src.at({ParameterType::BSM, "YL", 22})->get_val();

        dep_param->set_expected(-m_b * m_tau * ld * l_tau / std::pow(mH, 2));
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(4051516, 3231, 0, 1)}, sources, func);

    // double m_b = QCDHelper::mass_b_msbar();
    // double m_tau = sm("MASS", 15);
    // double l_tau = (*mod)("YL", 22);
    // // return this->double_to_complex_save("LO", -m_b * m_tau * thdm_params->ld * l_tau / std::pow(thdm_params->m_H, 2));
}

void C_S2_THDM::LO_calculation() {

    std::unordered_set<ParamId> sources {
        {"WPARAM_SI_SM", 4},
        {"WPARAM_SI_BSM", 7},
        {"WPARAM_SI_BSM", 8},
        {ParameterType::SM, "MASS", 15},
        {ParameterType::BSM, "MASS", 37},
        {ParameterType::BSM, "YL", {2, 2}},
        {ParameterType::SM, "MASS", 4},
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
        double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
        double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
        double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();

        double m_c = src.at({ParameterType::SM, "MASS", 4})->get_val(); // TODO : mass c does not run ?
        double m_tau = src.at({ParameterType::SM, "MASS", 15})->get_val();
        double l_tau = src.at({ParameterType::BSM, "YL", {2, 2}})->get_val();

        dep_param->set_expected(-m_c * m_tau * lu * l_tau / std::pow(mH, 2));
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(4051516, 3131, 0, 1)}, sources, func);

    // double m_c = sm("MASS", 4);
    // double m_tau = sm("MASS", 15);
    // double l_tau = (*mod)("YL", 22);
    // // return this->double_to_complex_save("LO", -m_c * m_tau * thdm_params->lu * l_tau / std::pow(thdm_params->m_H, 2));
}

void WilsonCoefficient_THDM::init(QCDOrder order) {
    if (!is_owned) {
        thdm_parameters::init();
    }

    WilsonCoefficient::init(order);
}

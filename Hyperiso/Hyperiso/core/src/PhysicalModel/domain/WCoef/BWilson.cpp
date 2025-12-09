#include "BWilson.h"

// ---------- C1 ----------

C1::C1() : WilsonCoefficient("C1", GroupMapper::str(WGroup::B, ScaleType::MATCHING)) {

    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
    matching_info[QCDOrder::NLO] = {
        {{"WPARAM_MATCH_SM", 3}},               // L
        [](const auto& src) {
            auto L = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 3);
            return 15. + 6. * L;
        },
        get_lhaid_from_name(QCDOrder::NLO)
    };

    matching_info[QCDOrder::NNLO] = {
        {
            {"WPARAM_MATCH_SM", 3},               // L
            {"WPARAM_MATCH_SM", LhaID(2, 1)}     // xt
        },
        [](const auto& src) {
            auto L = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 3);
            auto xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
            return -T(xt) + 7987./72. + 17. * PI2 / 3. + 475./6. * L + 17. * L * L;
        },
        get_lhaid_from_name(QCDOrder::NNLO)
    };
}

// ---------- C2 ----------

C2::C2() : WilsonCoefficient("C2", GroupMapper::str(WGroup::B, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {}, 
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };

    matching_info[QCDOrder::NLO] = MatchingInfo(LhaID(3040405, 4141, 1, 0));

    matching_info[QCDOrder::NNLO] = {
        {
            {"WPARAM_MATCH_SM", 3}               // L
        },
        compute_NNLO,
        get_lhaid_from_name(QCDOrder::NNLO)
    };
}

double C2::compute_LO(const ParamSrc& src) {
    return 1.;
}

double C2::compute_NNLO(const ParamSrc& src) {
    auto L = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 3);
    return 127. / 18. + 4. / 3. * PI2 + 46. / 3. * L + 4. * L * L;
}


// ---------- C3 ----------

C3::C3() : WilsonCoefficient("C3", GroupMapper::str(WGroup::B, ScaleType::MATCHING)) {

    matching_info[QCDOrder::NLO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));

    matching_info[QCDOrder::NLO] = MatchingInfo(get_lhaid_from_name(QCDOrder::NLO));

    matching_info[QCDOrder::NNLO] = {
        {
            {"WPARAM_MATCH_SM", 3},               // L
            {"WPARAM_MATCH_SM", 6},               // mass_top_muW
            {"WPARAM_MATCH_SM", LhaID(2, 1)},     // xt
            {"EW_SCALE", 1}                       // Q_match
        },
        compute_NNLO,
        get_lhaid_from_name(QCDOrder::NNLO)
    };
}

double C3::compute_NNLO(const ParamSrc& src) {
    double L = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 3);
    double mass_top_muW = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 6);
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    double Q_match = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);

    double coeff_temp = G1t(xt, log(Q_match * Q_match / (mass_top_muW * mass_top_muW)))
                      - 680. / 243.
                      - 20. / 81. * PI2
                      - 68. / 81. * L
                      - 20. / 27. * L * L;

    return coeff_temp;
}

// ---------- C4 ----------

C4::C4() : WilsonCoefficient("C4", GroupMapper::str(WGroup::B, ScaleType::MATCHING)) {

    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));

    // NLO
    matching_info[QCDOrder::NLO] = {
        {
            {"WPARAM_MATCH_SM", 3},               // L
            {"WPARAM_MATCH_SM", LhaID(2, 1)}      // x_t
        },
        compute_NLO,
        get_lhaid_from_name(QCDOrder::NLO)
    };

    // NNLO
    matching_info[QCDOrder::NNLO] = {
        {
            {"WPARAM_MATCH_SM", 3},               // L
            {"WPARAM_MATCH_SM", 6},               // mass_top_muW
            {"WPARAM_MATCH_SM", LhaID(2, 1)},     // x_t
            {"EW_SCALE", 1}                       // Q_match
        },
        compute_NNLO,
        get_lhaid_from_name(QCDOrder::NNLO)
    };
}

double C4::compute_NLO(const ParamSrc& src) {
    double L = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 3);
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    return E0t(xt) - 7. / 9. + 2. / 3. * L;
}

double C4::compute_NNLO(const ParamSrc& src) {
    double L = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 3);
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    double mass_top_muW = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 6);
    double Q_match = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);
    return E1t(xt, log(Q_match * Q_match / (mass_top_muW * mass_top_muW)))
           + 950. / 243.
           + 10. / 81. * PI2
           + 124. / 27. * L
           + 10. / 27. * L * L;
}

// ---------- C5 ----------

C5::C5() : WilsonCoefficient("C5", GroupMapper::str(WGroup::B, ScaleType::MATCHING)) {

    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));

    matching_info[QCDOrder::NLO] = MatchingInfo(get_lhaid_from_name(QCDOrder::NLO));

    matching_info[QCDOrder::NNLO] = {
        {
            {"WPARAM_MATCH_SM", 3},               // L
            {"WPARAM_MATCH_SM", 6},               // mass_top_muW
            {"WPARAM_MATCH_SM", LhaID(2, 1)},     // xt
            {"EW_SCALE", 1}                       // Q_match
        },
        compute_NNLO,
        get_lhaid_from_name(QCDOrder::NNLO)
    };
}

double C5::compute_NNLO(const ParamSrc& src) {
    double L = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 3);
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    double mass_top_muW = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 6);
    double Q_match = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);

    double coeff_temp =
        -G1t(xt, log(Q_match * Q_match / (mass_top_muW * mass_top_muW))) / 10.
        + 2. / 15. * E0t(xt)
        + 68. / 243.
        + 2. / 81. * PI2
        + 14. / 81. * L
        + 2. / 27. * L * L;

    return coeff_temp;
}

// ---------- C6 ----------

C6::C6() : WilsonCoefficient("C6", GroupMapper::str(WGroup::B, ScaleType::MATCHING)) {

    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));

    matching_info[QCDOrder::NLO] = MatchingInfo(get_lhaid_from_name(QCDOrder::NLO));

    matching_info[QCDOrder::NNLO] = {
        {
            {"WPARAM_MATCH_SM", 3},               // L
            {"WPARAM_MATCH_SM", 6},               // mass_top_muW
            {"WPARAM_MATCH_SM", LhaID(2, 1)},     // xt
            {"EW_SCALE", 1}                       // Q_match
        },
        compute_NNLO,
        get_lhaid_from_name(QCDOrder::NNLO)
    };
}

double C6::compute_NNLO(const ParamSrc& src) {
    double L = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 3);
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    double mass_top_muW = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 6);
    double Q_match = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);

    double coeff_temp =
        -3. / 16. * G1t(xt, log(Q_match * Q_match / (mass_top_muW * mass_top_muW)))
        + E0t(xt) / 4.
        + 85. / 162.
        + 5. / 108. * PI2
        + 35. / 108. * L
        + 5. / 36. * L * L;

    return coeff_temp;
}


// ---------- C7 ----------

C7::C7() : WilsonCoefficient("C7", GroupMapper::str(WGroup::B, ScaleType::MATCHING)) {
    // LO
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2, 1)}  // x_t
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };

    // std::cout << "matching_info LO : " << get_lhaid_from_name(QCDOrder::LO) << std::endl;
    // NLO
    matching_info[QCDOrder::NLO] = {
        {
            {"WPARAM_MATCH_SM", 3},           // L
            {"WPARAM_MATCH_SM", 6},           // mass_top_muW
            {"WPARAM_MATCH_SM", LhaID(2, 1)}, // x_t
            {"EW_SCALE", 1}                   // Q_match
        },
        compute_NLO,
        get_lhaid_from_name(QCDOrder::NLO)
    };

    // NNLO
    matching_info[QCDOrder::NNLO] = {
        {
            {"WPARAM_MATCH_SM", 3},           // L
            {"WPARAM_MATCH_SM", 6},           // mass_top_muW
            {"WPARAM_MATCH_SM", LhaID(2, 1)}, // x_t
            {"WPARAM_MATCH_SM", 7},           // xtW (top pole)
            {"WPARAM_MATCH_SM", 8},           // xtt
            {"EW_SCALE", 1},                  // Q_match
            {ParameterType::SM, "MASS", 24}   // m_W
        },
        compute_NNLO,
        get_lhaid_from_name(QCDOrder::NNLO)
    };
}

double C7::compute_LO(const ParamSrc& src) {
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    return -0.5 * A0t(xt) - 23. / 36.;
}

double C7::compute_NLO(const ParamSrc& src) {
    double L = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 3);
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    double mass_top_muW = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 6);
    double Q_match = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);

    return -0.5 * A1t(xt, log(Q_match * Q_match / (mass_top_muW * mass_top_muW)))
           + 713. / 243. + 4. / 81. * L
           - 4. / 9. * (E0t(xt) - 7. / 9. + 2. / 3. * L);
}

double C7::compute_NNLO(const ParamSrc& src) {
    double L     = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 3);
    double xt    = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    double xtW   = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 7);
    double xtt   = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 8);
    double mtop  = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 6);
    double Q     = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);
    double mW    = src.get_val(ParameterType::SM, "MASS", 24);

    double logqt = log(Q * Q / (mtop * mtop));
    double logqW = log(Q * Q / (mW * mW));

    double coeff_temp = C7t2mt(xtt)
        + logqt * (
            ((-592. * pow(xt, 5.) - 22. * pow(xt, 4.) + 12814. * pow(xt, 3.) - 6376. * xt * xt + 512. * xt)
                / 27. / pow(xt - 1., 5.)) * Li2(1. - 1. / xt)
            + ((-26838. * pow(xt, 5.) + 25938. * pow(xt, 4.) + 627367. * pow(xt, 3.) - 331956. * xt * xt + 16989. * xt - 460.)
                / 729. / pow(xt - 1., 6.)) * log(xt)
            + ((34400. * pow(xt, 5.) + 276644. * pow(xt, 4.) - 2668324. * pow(xt, 3.) + 1694437. * xt * xt - 323354. * xt + 53077.)
                / 2187. / pow(xt - 1., 5.))
            + logqt * (
                ((-63. * pow(xt, 5.) + 532. * pow(xt, 4.) + 2089. * pow(xt, 3.) - 1118. * xt * xt)
                    / 9. / pow(xt - 1., 6.)) * log(xt)
                + ((1186. * pow(xt, 5.) - 2705. * pow(xt, 4.) - 24791. * pow(xt, 3.) - 16099. * xt * xt + 19229. * xt - 2740.)
                    / 162. / pow(xt - 1., 5.))
            )
        )
        - (C7c2MW(xtW) + 13763. / 2187. * logqW + 814. / 729. * pow(logqW, 2.));

    return coeff_temp;
}

// ---------- C8 ----------

C8::C8() : WilsonCoefficient("C8", GroupMapper::str(WGroup::B, ScaleType::MATCHING)) {
    // LO
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)}  // x_t
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };

    // NLO
    matching_info[QCDOrder::NLO] = {
        {
            {"WPARAM_MATCH_SM", 3},           // L
            {"WPARAM_MATCH_SM", 6},           // mass_top_muW
            {"WPARAM_MATCH_SM", LhaID(2, 1)}, // x_t
            {"EW_SCALE", 1}                   // Q_match
        },
        compute_NLO,
        get_lhaid_from_name(QCDOrder::NLO)
    };

    // NNLO
    matching_info[QCDOrder::NNLO] = {
        {
            {"WPARAM_MATCH_SM", 3},           // L
            {"WPARAM_MATCH_SM", 6},           // mass_top_muW
            {"WPARAM_MATCH_SM", LhaID(2, 1)}, // x_t
            {"WPARAM_MATCH_SM", 7},           // xtW
            {"WPARAM_MATCH_SM", 8},           // xtt
            {"EW_SCALE", 1},                  // Q_match
            {ParameterType::SM, "MASS", 24}   // m_W
        },
        compute_NNLO,
        get_lhaid_from_name(QCDOrder::NNLO)
    };
}

double C8::compute_LO(const ParamSrc& src) {
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    return -0.5 * F0t(xt) - 1. / 3.;
}

double C8::compute_NLO(const ParamSrc& src) {
    double L = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 3);
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    double mass_top_muW = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 6);
    double Q_match = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);

    return -0.5 * F1t(xt, log(Q_match * Q_match / (mass_top_muW * mass_top_muW)))
           + 91. / 324.
           - 4. / 27. * L
           - (E0t(xt) - 7. / 9. + 2. / 3. * L) / 6.;
}

double C8::compute_NNLO(const ParamSrc& src) {
    double L     = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 3);
    double xt    = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    double xtW   = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 7);
    double xtt   = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 8);
    double mtop  = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 6);
    double Q     = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);
    double mW    = src.get_val(ParameterType::SM, "MASS", 24);

    double logqt = log(Q * Q / (mtop * mtop));
    double logqW = log(Q * Q / (mW * mW));

    double coeff_temp = C8t2mt(xtt)
        + logqt * (
            ((-148. * pow(xt, 5.) + 1052. * pow(xt, 4.) - 4811. * pow(xt, 3.) - 3520. * xt * xt - 61. * xt)
                / 18. / pow(xt - 1., 5.)) * Li2(1. - 1. / xt)
            + ((-15984. * pow(xt, 5.) + 152379. * pow(xt, 4.) - 1358060. * pow(xt, 3.) - 1201653. * xt * xt - 74190. * xt + 9188.)
                / 1944. / pow(xt - 1., 6.)) * log(xt)
            + ((109669. * pow(xt, 5.) - 1112675. * pow(xt, 4.) + 6239377. * pow(xt, 3.) + 8967623. * xt * xt + 768722. * xt - 42796.)
                / 11664. / pow(xt - 1., 5.))
            + logqt * (
                ((-139. * pow(xt, 4.) - 2938. * pow(xt, 3.) - 2683. * xt * xt)
                    / 12. / pow(xt - 1., 6.)) * log(xt)
                + ((1295. * pow(xt, 5.) - 7009. * pow(xt, 4.) + 29495. * pow(xt, 3.) + 64513. * xt * xt + 17458. * xt - 2072.)
                    / 216. / pow(xt - 1., 5.))
            )
        )
        - (C8c2MW(xtW) + 16607. / 5832. * logqW + 397. / 486. * pow(logqW, 2.));

    return coeff_temp;
}

// ---------- C9 ----------


C9::C9() : WilsonCoefficient("C9", GroupMapper::str(WGroup::B, ScaleType::MATCHING)) {
    // LO
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)},  // x_t
            {"WPARAM_MATCH_SM", 3},           // L
            {"WPARAM_SI_SM", 4}               // sw2
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };

    // NLO
    matching_info[QCDOrder::NLO] = {
        {
            {"WPARAM_MATCH_SM", 3},           // L
            {"WPARAM_MATCH_SM", 6},           // mass_top_muW
            {"WPARAM_MATCH_SM", LhaID(2, 1)}, // x_t
            {"WPARAM_SI_SM", 4},              // sw2
            {"EW_SCALE", 1}                   // Q_match
        },
        compute_NLO,
        get_lhaid_from_name(QCDOrder::NLO)
    };

    matching_info[QCDOrder::NNLO] = MatchingInfo(get_lhaid_from_name(QCDOrder::NNLO));

}

double C9::compute_LO(const ParamSrc& src) {
    double L    = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 3);
    double xt   = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    double sw2  = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);

    // printf("L = %.5f\n", L);
    // printf("x_t = %.5f\n", xt);
    // printf("sw2 = %.5f\n", sw2);

    // printf("C9_LO = %.5f\n", (1. - 4. * sw2) / sw2 * C0t(xt)
    //      - B0t(xt) / sw2
    //      - D0t(xt)
    //      + 1. / (4. * sw2)
    //      + 38. / 27.
    //      - 4. / 9. * L);

    return (1. - 4. * sw2) / sw2 * C0t(xt)
         - B0t(xt) / sw2
         - D0t(xt)
         + 1. / (4. * sw2)
         + 38. / 27.
         - 4. / 9. * L;
}

double C9::compute_NLO(const ParamSrc& src) {
    double L    = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 3);
    double xt   = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    double mtop = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 6);
    double Q    = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);
    double sw2  = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);

    double logqt = log(Q * Q / (mtop * mtop));

    // printf("L = %.5f\n", L);
    // printf("x_t = %.5f\n", xt);
    // printf("mtop = %.5f\n", mtop);
    // printf("Q = %.5f\n", Q);
    // printf("sw2 = %.5f\n", sw2);
    // printf("logqt = %.5f\n", logqt);

    // printf("C9_NLO = %.5f\n", (1. - 4. * sw2) / sw2 * C1t(xt, logqt)
    //      - B1t(xt, logqt) / sw2
    //      - D1t(xt, logqt)
    //      + 1. / sw2
    //      + 524. / 729.
    //      - 128. / 243. * PI2
    //      - 16. / 3. * L
    //      - 128. / 81. * L * L);

    return (1. - 4. * sw2) / sw2 * C1t(xt, logqt)
         - B1t(xt, logqt) / sw2
         - D1t(xt, logqt)
         + 1. / sw2
         + 524. / 729.
         - 128. / 243. * PI2
         - 16. / 3. * L
         - 128. / 81. * L * L;
}

// ---------- C10 ----------

C10::C10() : WilsonCoefficient("C10", GroupMapper::str(WGroup::B, ScaleType::MATCHING)) {
    // LO
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)},  // x_t
            {"WPARAM_SI_SM", 4}               // sw2
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };

    // NLO
    matching_info[QCDOrder::NLO] = {
        {
            {"WPARAM_MATCH_SM", 6},           // mtop at muW
            {"WPARAM_MATCH_SM", LhaID(2, 1)}, // x_t
            {"WPARAM_SI_SM", 4},              // sw2
            {"EW_SCALE", 1}                   // Q_match
        },
        compute_NLO,
        get_lhaid_from_name(QCDOrder::NLO)
    };

    // NNLO
    matching_info[QCDOrder::NNLO] = {
        {
            {"WPARAM_MATCH_SM", 6},           // mtop_muW
            {"WPARAM_MATCH_SM", LhaID(2, 1)}, // x_t
            {"WPARAM_MATCH_SM", 7},           // xtW
            {"WPARAM_MATCH_SM", 8},           // xtt
            {"WPARAM_SI_SM", 4},              // sw2
            {"EW_SCALE", 1},                  // Q_match
            {ParameterType::SM, "MASS", 24}   // m_W
        },
        compute_NNLO,
        get_lhaid_from_name(QCDOrder::NNLO)
    };
}

double C10::compute_LO(const ParamSrc& src) {
    double xt  = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    double sw2 = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);

    // printf("x_t = %.5f\n", xt);
    // printf("sw2 = %.5f\n", sw2);

    return (B0t(xt) - C0t(xt) - 0.25) / sw2;
}

double C10::compute_NLO(const ParamSrc& src) {
    double xt    = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    double mtop  = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 6);
    double Q     = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);
    double sw2   = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);

    // printf("x_t = %.5f\n", xt);
    // printf("mtop = %.5f\n", mtop);
    // printf("mu_W = %.5f\n", Q);
    // printf("sw2 = %.5f\n", sw2);

    double logqt = log(Q * Q / (mtop * mtop));
    return (B1t(xt, logqt) - C1t(xt, logqt)) / sw2 - 1. / sw2;
}

double C10::compute_NNLO(const ParamSrc& src) {
    double xt    = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});
    double xtW   = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 7);
    double xtt   = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 8);
    double mtop  = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 6);
    double Q     = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);
    double mW    = src.get_val(ParameterType::SM, "MASS", 24);
    double sw2   = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);

    // printf("x_t = %.5f\n", xt);
    // printf("mtop = %.5f\n", mtop);
    // printf("mu_W = %.5f\n", Q);
    // printf("sw2 = %.5f\n", sw2);
    // printf("xtW = %.5f\n", xtW);
    // printf("xtt = %.5f\n", xtt);
    // printf("mW = %.5f\n", mW);

    double logqt = log(Q * Q / (mtop * mtop));
    double logqW = log(Q * Q / (mW * mW));

    double coeff_temp =
        (C10Wt2mt(xtt)
         + logqt * (
             (69. + 1292. * xt - 209. * xt * xt) / 18. / pow(xt - 1., 3.)
             - (521. * xt + 105. * xt * xt - 50. * pow(xt, 3.)) / 9. / pow(xt - 1., 4.) * log(xt)
             - (47. * xt + xt * xt) / 3. / pow(xt - 1., 3.) * Li2(1. - 1. / xt)
             + logqt * (
                 (61. * xt + 11. * xt * xt) / 3. / pow(xt - 1., 3.)
                 - (49. * xt + 96. * xt * xt - pow(xt, 3.)) / 6. / pow(xt - 1., 4.) * log(xt)
             )
         )
         - (C10Wc2MW(xtW) - 23. / 6. * logqW)
         + C10Zt2mt(xtt)
         + logqt * (
             (188. * xt + 4. * xt * xt + 95. * pow(xt, 3.) - 47. * pow(xt, 4.)) / 6. / pow(xt - 1., 3.) * Li2(1. - 1. / xt)
             + (1468. * xt + 1578. * xt * xt - 25. * pow(xt, 3.) - 141. * pow(xt, 4.)) / 18. / pow(xt - 1., 4.) * log(xt)
             - (4622. * xt + 1031. * xt * xt + 582. * pow(xt, 3.) - 475. * pow(xt, 4.)) / 36. / pow(xt - 1., 3.)
             + logqt * (
                 (49. * xt + 315. * xt * xt - 4. * pow(xt, 3.)) / 6. / pow(xt - 1., 4.) * log(xt)
                 - (440. * xt + 257. * xt * xt + 72. * pow(xt, 3.) - 49. * pow(xt, 4.)) / 12. / pow(xt - 1., 3.)
             )
         )
         + C10Z2tri(xtt)
        ) * (-2. / sw2);

    return coeff_temp;
}



#include "BWilsonSUSY.h"
#include "NMSSMScalarMatching.h"

C1_susy::C1_susy() : WilsonCoefficient("C1_SUSY", GroupMapper::str(WGroup::B, ScaleType::MATCHING)) {
    matching_info[QCDOrder::NNLO] = {
        {
            {ParameterType::SM,    "MASS", 24},         // mW
            {ParameterType::BSM,   "MASS", 1000002},    // msq1
            {ParameterType::BSM,   "MASS", 1000004},    // msq2
            {ParameterType::BSM,   "MASS", 1000006},    // msq3
            {ParameterType::BSM,   "MASS", 2000002},    // msq4
            {ParameterType::BSM,   "MASS", 2000004},    // msq5
            {ParameterType::BSM,   "MASS", 2000006},    // msq6
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 2}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 3}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 4}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 5}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {15, 0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {15, 1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {15, 2}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {15, 3}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {15, 4}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {15, 5}}
        },
        compute_NNLO,
        get_lhaid_from_name(QCDOrder::NNLO)
    };
}

scalar_t C1_susy::compute_NNLO(const ParamSrc& src) {
    complex_t C1squark_2 = 0.0;
    double mW = src.get_val(ParameterType::SM, "MASS", 24);

    Array1D_7 MsqU = {
        src.get_val(ParameterType::BSM, "MASS", 1000002),
        src.get_val(ParameterType::BSM, "MASS", 1000004),
        src.get_val(ParameterType::BSM, "MASS", 1000006),
        src.get_val(ParameterType::BSM, "MASS", 2000002),
        src.get_val(ParameterType::BSM, "MASS", 2000004),
        src.get_val(ParameterType::BSM, "MASS", 2000006)
    };

    if (std::all_of(begin(MsqU), end(MsqU), [&](double m) { return std::abs(m) > mW / 2.0; })) {
        C1squark_2 = -208.0 / 3.0;
        for (int ae = 0; ae < 6; ++ae) {
            for (int row : {14, 15}) {
                double xsqa = pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {row, ae}) / mW, 2.0);
                double angle = 2.0 * asin(0.5 / sqrt(xsqa));
                
                C1squark_2 += -2.0 * pow(4.0 * xsqa - 1.0, 1.5) * Cl2(angle);
                C1squark_2 += 8.0 * (xsqa - 1.0 / 3.0) * log(xsqa) + 16.0 * xsqa;
            }
        }
    }

    return C1squark_2;
}


C3_susy::C3_susy() : WilsonCoefficient("C3_SUSY", GroupMapper::str(WGroup::B, ScaleType::MATCHING)) {
    matching_info[QCDOrder::NNLO] = {
        {
            {ParameterType::WILSON, "EW_SCALE", 1},
            {ParameterType::SM, "MASS", 24},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13, 0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13, 1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 2}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 3}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 4}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 5}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 0, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 0, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 1, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 1, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 2, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 2, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 3, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 3, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 4, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 4, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 5, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 5, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 0, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 0, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 1, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 1, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 2, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 2, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 3, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 3, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 4, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 4, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 5, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 5, 2}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", 19},           // kappa
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},       // yt
            {ParameterType::WILSON, "WPARAM_SI_BSM", 7},           // lu
            {ParameterType::BSM, "MASS", 37}                       // mH
        },
        compute_NNLO,
        get_lhaid_from_name(QCDOrder::NNLO)
    };
}

scalar_t C3_susy::compute_NNLO(const ParamSrc& src) {
    complex_t C3charg_2 = 0.0;
    double Q_match = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);;
    double mW = src.get_val(ParameterType::SM, "MASS", 24);

    for (int ie = 0; ie < 2; ++ie) {
        double Mch = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie});
        double ratio_mass_W_Mch = std::pow(mW / Mch, 2.0);

        for (int ae = 0; ae < 6; ++ae) {
            double MsqU = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae});
            double ratio_MsqU_Mch = std::pow(MsqU / Mch, 2.0);
            double log_mu_W_MsqU = std::log(std::pow(Q_match / MsqU, 2.0));

            double aij = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, ie, ae, 1});
            double bij = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, ie, ae, 2});

            C3charg_2 += ratio_mass_W_Mch * aij * bij * h71(ratio_MsqU_Mch, log_mu_W_MsqU);
        }
    }

    C3charg_2 *= src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 19); // kappa

    double yt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    double lu = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 7);
    double mH = src.get_val(ParameterType::BSM, "MASS", 37);

    double C3H_2 = G3H(yt, lu) + Delta3H(yt, lu) * std::log(std::pow(Q_match / mH, 2.0));

    return scalar_t(C3charg_2 + C3H_2);
}

C4_susy::C4_susy() : WilsonCoefficient("C4_SUSY", GroupMapper::str(WGroup::B, ScaleType::MATCHING)) {
    // NLO Matching
    matching_info[QCDOrder::NLO] = {
        {
            {ParameterType::WILSON, "WPARAM_SI_BSM", 7},        // lu
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},    // yt
            {ParameterType::WILSON, "EW_SCALE", 1},            // Q_match
            {ParameterType::BSM, "MASS", 37},                  // mH
            {ParameterType::SM, "MASS", 24},                   // mW
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13, 0}}, // Mch
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13, 1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 0}}, // MsqU
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 2}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 3}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 4}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 5}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 0, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 0, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 1, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 1, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 2, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 2, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 3, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 3, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 4, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 4, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 5, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 5, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 0, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 0, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 1, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 1, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 2, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 2, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 3, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 3, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 4, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 4, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 5, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 5, 2}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", 19}         // kappa
        },
        compute_NLO,
        get_lhaid_from_name(QCDOrder::NLO)
    };

    auto nnlo_sources = matching_info[QCDOrder::NLO].sources;
    for (int ae = 0; ae < 6; ++ae) {
        for (int be = 0; be < 6; ++be) {
            nnlo_sources.insert({ParameterType::WILSON, "MATRIX_BSM", {9, ae, be}});
        }
    }

    // NNLO Matching
    matching_info[QCDOrder::NNLO] = {
        nnlo_sources,
        compute_NNLO,
        get_lhaid_from_name(QCDOrder::NNLO)
    };
}

scalar_t C4_susy::compute_NLO(const ParamSrc& src) {
    complex_t C4charg_1 = 0.0;
    double mW = src.get_val(ParameterType::SM, "MASS", 24);

    for (int ie = 0; ie < 2; ++ie) {
        double Mch = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie});
        double ratio = std::pow(mW / Mch, 2.0);

        for (int ae = 0; ae < 6; ++ae) {
            double MsqU = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae});
            double x = std::pow(MsqU / Mch, 2.0);
            double aij = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, ie, ae, 1});
            double bij = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, ie, ae, 2});

            C4charg_1 += ratio * aij * bij * h40(x);
        }
    }

    C4charg_1 *= src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 19); // kappa

    double yt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    double lu = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 7);
    complex_t C4H_1 = EH(yt, lu);

    return scalar_t(C4H_1 + C4charg_1);
}

scalar_t C4_susy::compute_NNLO(const ParamSrc& src) {
    complex_t C4charg_2 = 0.0;
    complex_t C4four_2 = 0.0;

    double lu = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 7);
    double yt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    double Q_match = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);;
    double mW = src.get_val(ParameterType::SM, "MASS", 24);
    double mH = src.get_val(ParameterType::BSM, "MASS", 37);

    for (int ie = 0; ie < 2; ++ie) {
        double Mch = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie});
        double ratio_W_Mch = std::pow(mW / Mch, 2.0);

        for (int ae = 0; ae < 6; ++ae) {
            double MsqU_ae = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae});
            double x = std::pow(MsqU_ae / Mch, 2.0);
            double log_mu = std::log(std::pow(Q_match / MsqU_ae, 2.0));

            double aij = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, ie, ae, 1});
            double bij = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, ie, ae, 2});

            C4charg_2 += ratio_W_Mch * aij * bij * h41(x, log_mu);

            for (int be = 0; be < 6; ++be) {
                double MsqU_be = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, be});
                double M9_ae_be = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9, ae, be});
                double M9_be_ce = 0.0;

                for (int ce = 0; ce < 6; ++ce) {
                    double MsqU_ce = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce});
                    double y = std::pow(MsqU_ce / Mch, 2.0);
                    M9_be_ce = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9, be, ce});
                    double aij_ce = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, ie, ce, 2});

                    C4four_2 += ratio_W_Mch * M9_ae_be * (MsqU_be / Mch) * M9_be_ce *
                                (1.0 + log_mu) * aij * aij_ce * q61(x, y);
                }
            }
        }
    }

    C4charg_2 *= src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 19);
    C4four_2 *= src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 19);

    double C4H_2 = G4H(yt, lu) + Delta4H(yt, lu) * std::log(std::pow(Q_match / mH, 2.0));

    return scalar_t(C4charg_2 + C4four_2 + C4H_2);
}

C5_susy::C5_susy() : WilsonCoefficient("C5_SUSY", GroupMapper::str(WGroup::B, ScaleType::MATCHING)) {
    matching_info[QCDOrder::NNLO] = {
        {
            {ParameterType::WILSON, "WPARAM_SI_BSM", 7},         // lu
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},     // yt
            {ParameterType::WILSON, "EW_SCALE", 1},             // Q_match
            {ParameterType::BSM, "MASS", 37},                   // mH
            {ParameterType::SM, "MASS", 24},                    // mW
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13, 0}},  // Mch
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13, 1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 0}},  // MsqU
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 2}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 3}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 4}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 5}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 0, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 0, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 1, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 1, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 2, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 2, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 3, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 3, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 4, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 4, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 5, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 5, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 0, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 0, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 1, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 1, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 2, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 2, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 3, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 3, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 4, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 4, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 5, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 5, 2}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", 19} 
        },
        compute_NNLO,
        get_lhaid_from_name(QCDOrder::NNLO)
    };
}

scalar_t C5_susy::compute_NNLO(const ParamSrc& src) {
    double lu      = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 7);
    double yt      = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    double Q_match = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);;
    double mW      = src.get_val(ParameterType::SM, "MASS", 24);
    double mH      = src.get_val(ParameterType::BSM, "MASS", 37);

    complex_t C3charg_2 = 0.0;
    complex_t C4charg_1 = 0.0;

    for (int ie = 0; ie < 2; ++ie) {
        double Mch = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie});

        for (int ae = 0; ae < 6; ++ae) {
            double MsqU = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae});
            double x    = std::pow(MsqU / Mch, 2.0);
            double logx = std::log(std::pow(Q_match / MsqU, 2.0));
            double aij  = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, ie, ae, 1});
            double bij  = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, ie, ae, 2});

            double ratio_mass_W_Mch = std::pow(mW / Mch, 2.0);

            C3charg_2 += ratio_mass_W_Mch * aij * bij * h71(x, logx);
            C4charg_1 += ratio_mass_W_Mch * aij * bij * h40(x);
        }
    }

    double kappa = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 19);
    C3charg_2 *= kappa;
    C4charg_1 *= kappa;

    complex_t C5charg_2 = -C3charg_2 / 10.0 + (2.0 / 15.0) * C4charg_1;

    double C4H_1 = EH(yt, lu);
    double C3H_2 = G3H(yt, lu) + Delta3H(yt, lu) * std::log(std::pow(Q_match / mH, 2.0));
    double C5H_2 = -C3H_2 / 10.0 + (2.0 / 15.0) * C4H_1;

    return scalar_t(C5charg_2 + C5H_2);
}

C6_susy::C6_susy() : WilsonCoefficient("C6_SUSY", GroupMapper::str(WGroup::B, ScaleType::MATCHING)) {
    matching_info[QCDOrder::NNLO] = {
        {
            {ParameterType::WILSON, "WPARAM_SI_BSM", 7},         // lu
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},     // yt
            {ParameterType::WILSON, "EW_SCALE", 1},             // Q_match
            {ParameterType::BSM, "MASS", 37},                   // mH
            {ParameterType::SM, "MASS", 24},                    // mW
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13, 0}},  // Mch
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13, 1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 0}},  // MsqU
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 2}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 3}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 4}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 5}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 0, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 0, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 1, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 1, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 2, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 2, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 3, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 3, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 4, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 4, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 5, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 5, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 0, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 0, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 1, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 1, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 2, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 2, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 3, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 3, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 4, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 4, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 5, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 5, 2}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", 19}          // kappa
        },
        compute_NNLO,
        get_lhaid_from_name(QCDOrder::NNLO)
    };
}

scalar_t C6_susy::compute_NNLO(const ParamSrc& src) {
    double lu      = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 7);
    double yt      = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    double Q_match = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);;
    double mW      = src.get_val(ParameterType::SM, "MASS", 24);
    double mH      = src.get_val(ParameterType::BSM, "MASS", 37);

    complex_t C3charg_2 = 0.0;
    complex_t C4charg_1 = 0.0;

    for (int ie = 0; ie < 2; ++ie) {
        double Mch = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie});

        for (int ae = 0; ae < 6; ++ae) {
            double MsqU = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae});
            double x    = std::pow(MsqU / Mch, 2.0);
            double logx = std::log(std::pow(Q_match / MsqU, 2.0));

            double aij = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, ie, ae, 1});
            double bij = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, ie, ae, 2});

            double ratio_mass_W_Mch = std::pow(mW / Mch, 2.0);

            C3charg_2 += ratio_mass_W_Mch * aij * bij * h71(x, logx);
            C4charg_1 += ratio_mass_W_Mch * aij * bij * h40(x);
        }
    }

    double kappa = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 19);
    C3charg_2 *= kappa;
    C4charg_1 *= kappa;

    complex_t C6charg_2 = -3.0 / 16.0 * C3charg_2 + 1.0 / 4.0 * C4charg_1;

    complex_t C4H_1 = EH(yt, lu);
    complex_t C3H_2 = G3H(yt, lu) + Delta3H(yt, lu) * std::log(std::pow(Q_match / mH, 2.0));
    complex_t C6H_2 = -3.0 / 16.0 * C3H_2 + 1.0 / 4.0 * C4H_1;

    return scalar_t(C6charg_2 + C6H_2);
}

C7_susy::C7_susy() : WilsonCoefficient("C7_SUSY", GroupMapper::str(WGroup::B, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}},
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},
            {ParameterType::SM, "MASS", 24},       
            // {ParameterType::WILSON, "EW_SCALE", 1},                   // mW
            {ParameterType::WILSON, "WPARAM_SI_BSM", 19},             // kappa
            {ParameterType::WILSON, "WPARAM_SI_BSM", 7},             // lu
            {ParameterType::WILSON, "WPARAM_SI_BSM", 8},             // ld
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13, 0}},       // Mch
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13, 1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 0}},       // MsqU
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 2}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 3}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 4}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 5}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 0, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 0, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 1, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 1, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 2, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 2, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 3, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 3, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 4, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 4, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 5, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 5, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 0, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 0, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 1, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 1, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 2, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 2, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 3, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 3, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 4, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 4, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 5, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 5, 2}},
            
            {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 0, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 0, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 1, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 1, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 2, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 2, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 3, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 3, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 4, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 4, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 5, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 5, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 0, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 0, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 1, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 1, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 2, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 2, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 3, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 3, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 4, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 4, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 5, 1}}, {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 5, 2}},
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},          // yt
            {ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}},       // xt
            {ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}},       // mb(Q)
            {ParameterType::WILSON, "EPSILON_SUSY", {0,1}},          // eps0
            {ParameterType::WILSON, "EPSILON_SUSY", {0,2}},          // eps0'
            {ParameterType::WILSON, "EPSILON_SUSY", 1},              // eps1'
            {ParameterType::WILSON, "EPSILON_SUSY", 2},              // eps2
            {ParameterType::WILSON, "EPSILON_SUSY", 3},              // epsb
            {ParameterType::WILSON, "EPSILON_SUSY", 4},              // epsb'
            {ParameterType::WILSON, "WPARAM_SI_BSM", 18},            // kappafactor
            {ParameterType::BSM, "MINPAR", 3},                         // tanb
            {ParameterType::BSM, "ALPHA", {}},                        // alpha
            {ParameterType::BSM, "MASS", 25},                        // mh
            {ParameterType::BSM, "MASS", 35},                        // mH0
            {ParameterType::BSM, "MASS", 36},                        // mA0
            // {ParameterType::BSM, "MASS", 37},                        // mH+
            {ParameterType::BSM, "MASS", 45},                       // mH03    
            {ParameterType::BSM, "MASS", 46},                       // mA02
            {ParameterType::BSM, "NMHMIX", {0+1, 0+1}}, {ParameterType::BSM, "NMHMIX", {0+1, 1+1}},
            {ParameterType::BSM, "NMHMIX", {1+1, 0+1}}, {ParameterType::BSM, "NMHMIX", {1+1, 1+1}},
            {ParameterType::BSM, "NMHMIX", {2+1, 0+1}}, {ParameterType::BSM, "NMHMIX", {2+1, 1+1}},
            {ParameterType::BSM, "NMAMIX" , {0+1, 0+1}}, {ParameterType::BSM, "NMAMIX" , {0+1, 1+1}},
            {ParameterType::BSM, "NMAMIX" , {1+1, 0+1}}, {ParameterType::BSM, "NMAMIX" , {1+1, 1+1}}
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };

    matching_info[QCDOrder::NLO] = {
        {
            // {ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}},
            {ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", 19},   // kappa
            {ParameterType::WILSON, "WPARAM_SI_BSM", 7},   // lu
            {ParameterType::WILSON, "WPARAM_SI_BSM", 8},   // ld
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", 1}, // yt
            {ParameterType::WILSON, "EW_SCALE", 1},     // Q_match
            {ParameterType::SM, "MASS", 24},               // mW
            // {ParameterType::BSM, "MASS", 25},              // mh
            // {ParameterType::BSM, "MASS", 35},              // mH0
            {ParameterType::BSM, "MASS", 37},              // mH
        
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 0}},       // MsqU
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 2}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 3}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 4}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 5}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,3,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,4,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,5,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,5,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,3,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,4,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,5,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,5,2}},
        
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,5,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,5,2}},
        
            {ParameterType::WILSON, "MATRIX_BSM", {9,0,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,0,3}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,0,4}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,0,5}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,1,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,1,3}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,1,4}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,1,5}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,2,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,2,3}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,2,4}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,2,5}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,3,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,3,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,3,3}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,3,4}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,3,5}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,4,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,4,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,4,3}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,4,4}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,4,5}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,5,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,5,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,5,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,5,3}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,5,4}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,5,5}}
        },
        compute_NLO,
        get_lhaid_from_name(QCDOrder::NLO)
    };
    matching_info[QCDOrder::NNLO] = {
        {
            {ParameterType::WILSON, "WPARAM_SI_BSM", 7},     // lu
            {ParameterType::WILSON, "WPARAM_SI_BSM", 8},     // ld
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},  // yt
            {ParameterType::WILSON, "EW_SCALE", 1},          // Q_match
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6}    // mass_top_muW
        },
        compute_NNLO,
        get_lhaid_from_name(QCDOrder::NNLO)
    };
}

scalar_t C7_susy::compute_LO(const ParamSrc& src) {

    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1});
    double lu = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 7);
    double ld = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 8);
    double yt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    // double Q_match = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);;
    double mW = src.get_val(ParameterType::SM, "MASS", 24);
    // double mH = src.get_val(ParameterType::BSM, "MASS", 37);
    double tanb = src.get_val(ParameterType::BSM, "MINPAR", 3);
    double alpha = src.get_val(ParameterType::BSM, "ALPHA", {});
    double mH0 = src.get_val(ParameterType::BSM, "MASS", 35);
    double mh = src.get_val(ParameterType::BSM, "MASS", 25);
    double mA0 = src.get_val(ParameterType::BSM, "MASS", 36);
    double mass_b_muW = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1});

    // std::cout << "mass_b_muW : " << mass_b_muW << std::endl;
    // std::cout << "xt : " << xt << std::endl;
    // std::cout << "yt : " << yt << std::endl;
    // std::cout << "lu : " << lu << std::endl;
    // std::cout << "ld : " << ld << std::endl;

    double epsilon0 = src.get_val(ParameterType::WILSON, "EPSILON_SUSY", {0,1});
    double epsilon0p = src.get_val(ParameterType::WILSON, "EPSILON_SUSY", {0,2});
    double epsilon1p = src.get_val(ParameterType::WILSON, "EPSILON_SUSY", 1);
    double epsilon2 = src.get_val(ParameterType::WILSON, "EPSILON_SUSY", 2);
    double epsilonb = src.get_val(ParameterType::WILSON, "EPSILON_SUSY", 3);
    double epsilonbp = src.get_val(ParameterType::WILSON, "EPSILON_SUSY", 4);

    double m_H03 = src.get_val(ParameterType::BSM, "MASS", 45);
    double m_A02 = src.get_val(ParameterType::BSM, "MASS", 46);

    complex_t C7SMeps_0 = (epsilonb - epsilonbp) / (1. + epsilonb * tanb) * tanb * F7_2(xt);
    complex_t C7Heps_0 = (-epsilon0p - epsilonb) / (1. + epsilonb * tanb) * tanb * F7_2(yt);

    // std::cout << "epsilon0" << epsilon0 << std::endl;
    // std::cout << "epsilon2" << epsilon2 << std::endl;
    // std::cout << "epsilonb" << epsilonb << std::endl;
    // std::cout << "epsilon0p" << epsilon0p << std::endl;
    // std::cout << "epsilonbp" << epsilonbp << std::endl;
    // std::cout << "epsilon1p" << epsilon1p << std::endl;
    // std::cout << "tanb" << tanb << std::endl;

    // std::cout << "alpha" << alpha << std::endl;
    // std::cout << "mh" << mh << std::endl;
    // std::cout << "mH0" << mH0 << std::endl;

    // std::cout << "mW" << mW << std::endl;
    complex_t C7Heps2_0 = 0.;

    // MAJ : Use MODSEL block instead
    if ((m_A02 == 0.) && (m_H03 == 0.)) {
        C7Heps2_0 = -epsilon2 * epsilon1p * pow(tanb, 2.) / (1. + epsilonb * tanb) / (1. + epsilon0 * tanb) * F7_2(yt);
        C7Heps2_0 += epsilon2 / pow(1. + epsilonb * tanb, 2.) * (1. + pow(tanb, 2.)) / (1. + epsilon0 * tanb) / 72. *
            ((cos(alpha) + sin(alpha) * tanb) * (-sin(alpha) + epsilonb * cos(alpha)) * pow(mass_b_muW / mh, 2.) +
             (sin(alpha) - cos(alpha) * tanb) * (cos(alpha) + epsilonb * sin(alpha)) * pow(mass_b_muW / mH0, 2.) +
             (-cos(atan(tanb)) - sin(atan(tanb)) * tanb) * (sin(atan(tanb)) - epsilonb * cos(atan(tanb))) * pow(mass_b_muW / mA0, 2.));
    } else {
        C7Heps2_0 = -epsilon2 * epsilon1p * pow(tanb, 2.) / (1. + epsilonb * tanb) / (1. + epsilon0 * tanb) * F7_2(yt);
        C7Heps2_0 += epsilon2 / pow(1. + epsilonb * tanb, 2.) * (1. + pow(tanb, 2.)) / (1. + epsilon0 * tanb) / 72. *
            ((src.get_val(ParameterType::BSM, "NMHMIX", {0+1, 0+1}) + src.get_val(ParameterType::BSM, "NMHMIX", {0+1, 1+1}) * tanb) *
             (-src.get_val(ParameterType::BSM, "NMHMIX", {0+1, 1+1}) + epsilonb * src.get_val(ParameterType::BSM, "NMHMIX", {0+1, 0+1})) * pow(mass_b_muW / mh, 2.) +
             (src.get_val(ParameterType::BSM, "NMHMIX", {1+1, 0+1}) + src.get_val(ParameterType::BSM, "NMHMIX", {1+1, 1+1}) * tanb) *
             (-src.get_val(ParameterType::BSM, "NMHMIX", {1+1, 1+1}) + epsilonb * src.get_val(ParameterType::BSM, "NMHMIX", {1+1, 0+1})) * pow(mass_b_muW / mH0, 2.) +
             (src.get_val(ParameterType::BSM, "NMHMIX", {2+1, 0+1}) + src.get_val(ParameterType::BSM, "NMHMIX", {2+1, 1+1}) * tanb) *
             (-src.get_val(ParameterType::BSM, "NMHMIX", {2+1, 1+1}) + epsilonb * src.get_val(ParameterType::BSM, "NMHMIX", {2+1, 0+1})) * pow(mass_b_muW / m_H03, 2.) +
             (src.get_val(ParameterType::BSM, "NMAMIX", {0+1, 0+1}) + src.get_val(ParameterType::BSM, "NMAMIX", {0+1, 1+1}) * tanb) *
             (-src.get_val(ParameterType::BSM, "NMAMIX", {0+1, 1+1}) + epsilonb * src.get_val(ParameterType::BSM, "NMAMIX", {0+1, 0+1})) * pow(mass_b_muW / mA0, 2.) +
             (src.get_val(ParameterType::BSM, "NMAMIX", {1+1, 0+1}) + src.get_val(ParameterType::BSM, "NMAMIX", {1+1, 1+1}) * tanb) *
             (-src.get_val(ParameterType::BSM, "NMAMIX", {1+1, 1+1}) + epsilonb * src.get_val(ParameterType::BSM, "NMAMIX", {1+1, 0+1})) * pow(mass_b_muW / m_A02, 2.));
    }

    complex_t C7charg_0 = 0.0;
    complex_t C7_chargeps_0 = 0.0;

    auto calculateContribution = [src, mW, tanb, epsilonb](auto hFunc, int X, int X2, int ie, int ae, bool isChargeps) -> double {
        double Mch = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie});
        double ratio = std::pow(mW / Mch, 2);
        double x = pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}) / Mch, 2);
        double factor = isChargeps ? (-epsilonb / (1.0 + epsilonb * tanb) * tanb) : 1.0;

        return ratio *
            src.get_val(ParameterType::WILSON, "MATRIX_BSM", {X, ie, ae, 1}) *
            src.get_val(ParameterType::WILSON, "MATRIX_BSM", {X2, ie, ae, 2}) *
            hFunc(x) *
            src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 18) *
            factor;
    };
    // std::cout << "mass_b_muW : " << mass_b_muW << std::endl; 
    for (int ie = 0; ie < 2; ++ie) {
        for (int ae = 0; ae < 6; ++ae) {
            double Mch = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie});
            C7charg_0 += calculateContribution(h10, 3, 3, ie, ae, false);
            C7charg_0 += Mch / mass_b_muW * calculateContribution(h20, 3, 4, ie, ae, false);
            C7_chargeps_0 += Mch / mass_b_muW * calculateContribution(h20, 3, 4, ie, ae, true);
        }
    }

    complex_t C7H_0 = 1. / 3. * lu * lu * F7_1(yt) - lu * ld * F7_2(yt);

    return scalar_t(C7SMeps_0 + C7Heps_0 + C7Heps2_0 + C7charg_0 + C7_chargeps_0 + C7H_0);
}

scalar_t C7_susy::compute_NLO(const ParamSrc& src) {
    double lu      = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 7);
    double ld      = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 8);
    double yt      = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    double Q_match = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);;
    double mW      = src.get_val(ParameterType::SM, "MASS", 24);
    double mH      = src.get_val(ParameterType::BSM, "MASS", 37);
    double mass_b_muW = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1});

    complex_t C7charg_1 = 0.0;
    complex_t C7four_1 = 0.0;

    for (int ie = 0; ie < 2; ++ie) {
        double Mch_ie = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie});
        double ratio_W_Mch = std::pow(mW / Mch_ie, 2.0);

        for (int ae = 0; ae < 6; ++ae) {
            double MsqU_ae = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae});
            double x = std::pow(MsqU_ae / Mch_ie, 2.0);
            double log_mu = std::log(std::pow(Q_match / MsqU_ae, 2.0));

            double aij = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, ie, ae, 1});
            double bij = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, ie, ae, 2});
            double cij = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2});

            C7charg_1 += ratio_W_Mch * (aij * bij * h11(x, log_mu) + Mch_ie / mass_b_muW * aij * cij * h21(x, log_mu));

            for (int ce = 0; ce < 6; ++ce) {
                double MsqU_ce = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce});
                double log_ce = std::log(std::pow(Q_match / MsqU_ce, 2.0));
                double ratio_ce_ie = MsqU_ce / Mch_ie;

                for (int de = 0; de < 6; ++de) {
                    double MsqU_de = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de});
                    double x_ae = std::pow(MsqU_ae / Mch_ie, 2.0);
                    double x_de = std::pow(MsqU_de / Mch_ie, 2.0);

                    double Nij = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9, ae, ce});
                    double Ncd = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9, ce, de});

                    double bde = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, ie, de, 2});
                    double cde = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4, ie, de, 2});

                    C7four_1 += ratio_W_Mch * Nij * ratio_ce_ie * Ncd * (1.0 + log_ce) * (
                        aij * bde * (-q11(x_ae, x_de) + 2.0 / 3.0 * q21(x_ae, x_de)) +
                        Mch_ie / mass_b_muW * aij * cde * (-q31(x_ae, x_de) + 2.0 / 3.0 * q41(x_ae, x_de))
                    );
                }
            }
        }
    }

    C7charg_1 *= -0.5 * src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 19);  // kappa
    C7four_1 *= -0.5 * src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 19);  // kappa
    complex_t C7H_1 = G7H(yt, lu, ld) + Delta7H(yt, lu, ld) * std::log(std::pow(Q_match / mH, 2.0)) - 4.0 / 9.0 * EH(yt, lu);

    return scalar_t(C7charg_1 + C7four_1 + C7H_1);
}

scalar_t C7_susy::compute_NNLO(const ParamSrc& src) {
    double lu             = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 7);
    double ld             = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 8);
    double yt             = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    double Q_match        = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);;
    double mass_top_muW   = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 6);

    double log_term = std::log(std::pow(Q_match / mass_top_muW, 2.0));

    double coeff_temp = C7H2(yt, lu, ld, log_term);

    return scalar_t(coeff_temp);
}

C8_susy::C8_susy() : WilsonCoefficient("C8_SUSY", GroupMapper::str(WGroup::B, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_SI_BSM", 7},      // lu
            {ParameterType::WILSON, "WPARAM_SI_BSM", 8},      // ld
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},   // yt
            // {ParameterType::WILSON, "EW_SCALE", 1},           // Q_match
            {ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}}, // xt
            {ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}}, // mass_b_muW
            {ParameterType::SM, "MASS", 24},                  // mW
            {ParameterType::BSM, "MASS", 25},                 // mh
            {ParameterType::BSM, "MASS", 35},                 // mH0
            {ParameterType::BSM, "MASS", 36},                 // mA0
            // {ParameterType::BSM, "MASS", 37},                 // mH
            {ParameterType::BSM, "MASS", 45},                       // mH03    
            {ParameterType::BSM, "MASS", 46},                       // mA02
            {ParameterType::BSM, "MINPAR", 3},                  // tanb
            {ParameterType::BSM, "ALPHA", {}},                 // alpha
        
            // EPSILON_SUSY parameters
            {ParameterType::WILSON, "EPSILON_SUSY", {0,1}},
            {ParameterType::WILSON, "EPSILON_SUSY", {0,2}},
            {ParameterType::WILSON, "EPSILON_SUSY", 1},
            {ParameterType::WILSON, "EPSILON_SUSY", 2},
            {ParameterType::WILSON, "EPSILON_SUSY", 3},
            {ParameterType::WILSON, "EPSILON_SUSY", 4},
        
            // HMIX & AMIX matrix elements used in the else branch
            {ParameterType::BSM, "NMHMIX", {0+1, 0+1}},
            {ParameterType::BSM, "NMHMIX", {0+1, 1+1}},
            {ParameterType::BSM, "NMHMIX", {1+1, 0+1}},
            {ParameterType::BSM, "NMHMIX", {1+1, 1+1}},
            {ParameterType::BSM, "NMHMIX", {2+1, 0+1}},
            {ParameterType::BSM, "NMHMIX", {2+1, 1+1}},
            // {ParameterType::BSM, "NMHMIX", {3+1, 0+1}},
            {ParameterType::BSM, "NMAMIX" , {0+1, 0+1}},
            {ParameterType::BSM, "NMAMIX" , {0+1, 1+1}},
            {ParameterType::BSM, "NMAMIX" , {1+1, 0+1}},
            {ParameterType::BSM, "NMAMIX" , {1+1, 1+1}},
        
            // Additional WPARAM_SI_BSM entries for calculateContribution
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,2}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,3}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,4}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,5}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", 18},     // kappafactor
        
            // Matrices MATRIX_BSM (h50 and h60 contributions)
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,3,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,4,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,5,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,5,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,3,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,4,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,5,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,5,2}},
        
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,5,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,5,2}}
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };

    matching_info[QCDOrder::NLO] = {
        {
            {ParameterType::WILSON, "WPARAM_SI_BSM", 7},     // lu
            {ParameterType::WILSON, "WPARAM_SI_BSM", 8},     // ld
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},  // yt
            {ParameterType::WILSON, "EW_SCALE", 1},          // Q_match
            // {ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}}, // xt
            {ParameterType::SM, "MASS", 24},                 // mW
            {ParameterType::BSM, "MASS", 37},                // mH
            {ParameterType::BSM, "MINPAR", 3},                 // tanb
            {ParameterType::BSM, "ALPHA", {}},                // alpha
            {ParameterType::BSM, "MASS", 35},                // mH0
            {ParameterType::BSM, "MASS", 25},                // mh
            {ParameterType::BSM, "MASS", 36},                // mA0
            {ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}}, // mass_b_muW
            {ParameterType::WILSON, "WPARAM_SI_BSM", 19},     //kappa
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,2}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,3}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,4}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,5}},
        
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,3,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,4,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,5,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,5,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,3,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,4,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,5,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,5,2}},
        
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,5,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,5,2}},
        
            {ParameterType::WILSON, "MATRIX_BSM", {9,0,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,0,3}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,0,4}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,0,5}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,1,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,1,3}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,1,4}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,1,5}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,2,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,2,3}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,2,4}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,2,5}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,3,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,3,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,3,3}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,3,4}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,3,5}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,4,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,4,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,4,3}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,4,4}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,4,5}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,5,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,5,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,5,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,5,3}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,5,4}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,5,5}}
        },
        compute_NLO,
        get_lhaid_from_name(QCDOrder::NLO)
    };

    matching_info[QCDOrder::NNLO] = {
        {
            {ParameterType::WILSON, "WPARAM_SI_BSM", 7},     // lu
            {ParameterType::WILSON, "WPARAM_SI_BSM", 8},     // ld
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},  // yt
            {ParameterType::WILSON, "EW_SCALE", 1},          // Q_match
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6}    // mass_top_muW
        },
        compute_NNLO,
        get_lhaid_from_name(QCDOrder::NNLO)
    };

}

scalar_t C8_susy::compute_LO(const ParamSrc& src) {
    double xt      = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1});
    double lu      = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 7);
    double ld      = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 8);
    double yt      = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    // double Q_match = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);;
    double mW      = src.get_val(ParameterType::SM, "MASS", 24);
    // double mH      = src.get_val(ParameterType::BSM, "MASS", 37);
    double tanb    = src.get_val(ParameterType::BSM, "MINPAR", 3);
    double alpha   = src.get_val(ParameterType::BSM, "ALPHA", {});
    double mH0     = src.get_val(ParameterType::BSM, "MASS", 35);
    double mh      = src.get_val(ParameterType::BSM, "MASS", 25);
    double mA0     = src.get_val(ParameterType::BSM, "MASS", 36);
    double mb_muW  = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1});

    double eps0  = src.get_val(ParameterType::WILSON, "EPSILON_SUSY", {0,1});
    double eps0p = src.get_val(ParameterType::WILSON, "EPSILON_SUSY", {0,2});
    double eps1p = src.get_val(ParameterType::WILSON, "EPSILON_SUSY", 1);
    double eps2  = src.get_val(ParameterType::WILSON, "EPSILON_SUSY", 2);
    double epsb  = src.get_val(ParameterType::WILSON, "EPSILON_SUSY", 3);
    double epsbp = src.get_val(ParameterType::WILSON, "EPSILON_SUSY", 4);

    double m_H03 = src.get_val(ParameterType::BSM, "MASS", 45);
    double m_A02 = src.get_val(ParameterType::BSM, "MASS", 46);

    complex_t C8SMeps_0 = (epsb - epsbp) / (1. + epsb * tanb) * tanb * F8_2(xt);
    complex_t C8Heps_0  = (-eps0p - epsb) / (1. + epsb * tanb) * tanb * F8_2(yt);

    complex_t C8Heps2_0 = 0.0;

    if ((m_A02 == 0.) && (m_H03 == 0.)) {
        C8Heps2_0 = -eps2 * eps1p * pow(tanb, 2.) / ((1. + epsb * tanb) * (1. + eps0 * tanb)) * F8_2(yt);
        C8Heps2_0 += eps2 / pow(1. + epsb * tanb, 2.) * (1. + pow(tanb, 2.)) / (1. + eps0 * tanb) / 72. *
            ((cos(alpha) + sin(alpha) * tanb) * (-sin(alpha) + epsb * cos(alpha)) * pow(mb_muW / mh, 2.) +
             (sin(alpha) - cos(alpha) * tanb) * (cos(alpha) + epsb * sin(alpha)) * pow(mb_muW / mH0, 2.) +
             (-cos(atan(tanb)) - sin(atan(tanb)) * tanb) * (sin(atan(tanb)) - epsb * cos(atan(tanb))) * pow(mb_muW / mA0, 2.));
    } else {
        C8Heps2_0 = -eps2 * eps1p * pow(tanb, 2.) / ((1. + epsb * tanb) * (1. + eps0 * tanb)) * F8_2(yt);
        C8Heps2_0 += -3. * eps2 / pow(1. + epsb * tanb, 2.) * (1. + pow(tanb, 2.)) / (1. + eps0 * tanb) / 72. *
            ((src.get_val(ParameterType::BSM, "NMHMIX", {0+1, 0+1}) + src.get_val(ParameterType::BSM, "NMHMIX", {0+1, 1+1}) * tanb) *
             (-src.get_val(ParameterType::BSM, "NMHMIX", {0+1, 1+1}) + epsb * src.get_val(ParameterType::BSM, "NMHMIX", {0+1, 0+1})) * pow(mb_muW / mh, 2.) +
             (src.get_val(ParameterType::BSM, "NMHMIX", {1+1, 0+1}) + src.get_val(ParameterType::BSM, "NMHMIX", {1+1, 1+1}) * tanb) *
             (-src.get_val(ParameterType::BSM, "NMHMIX", {1+1, 1+1}) + epsb * src.get_val(ParameterType::BSM, "NMHMIX", {1+1, 0+1})) * pow(mb_muW / mH0, 2.) +
             (src.get_val(ParameterType::BSM, "NMHMIX", {2+1, 0+1}) + src.get_val(ParameterType::BSM, "NMHMIX", {2+1, 1+1}) * tanb) *
             (-src.get_val(ParameterType::BSM, "NMHMIX", {2+1, 1+1}) + epsb * src.get_val(ParameterType::BSM, "NMHMIX", {2+1, 0+1})) * pow(mb_muW / m_H03, 2.) +
             (src.get_val(ParameterType::BSM, "NMAMIX", {0+1, 0+1}) + src.get_val(ParameterType::BSM, "NMAMIX", {0+1, 1+1}) * tanb) *
             (-src.get_val(ParameterType::BSM, "NMAMIX", {0+1, 1+1}) + epsb * src.get_val(ParameterType::BSM, "NMAMIX", {0+1, 0+1})) * pow(mb_muW / mA0, 2.) +
             (src.get_val(ParameterType::BSM, "NMAMIX", {1+1, 0+1}) + src.get_val(ParameterType::BSM, "NMAMIX", {1+1, 1+1}) * tanb) *
             (-src.get_val(ParameterType::BSM, "NMAMIX", {1+1, 1+1}) + epsb * src.get_val(ParameterType::BSM, "NMAMIX", {1+1, 0+1})) * pow(mb_muW / m_A02, 2.));
    }

    auto calculateContribution = [src, mW, tanb, epsb](auto hFunc, int X, int X2, int ie, int ae, bool isChargeps) -> double {
        double Mch = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie});
        double x = pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}) / Mch, 2.0);
        double ratio = std::pow(mW / Mch, 2.0);
        double factor = isChargeps ? (-epsb / (1. + epsb * tanb) * tanb) : 1.0;
        double kappa = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 18);

        return ratio *
               src.get_val(ParameterType::WILSON, "MATRIX_BSM", {X, ie, ae, 1}) *
               src.get_val(ParameterType::WILSON, "MATRIX_BSM", {X2, ie, ae, 2}) *
               hFunc(x) * kappa * factor;
    };

    complex_t C8charg_0 = 0.0;
    complex_t C8_chargeps_0 = 0.0;

    for (int ie = 0; ie < 2; ++ie) {
        double Mch = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie});
        for (int ae = 0; ae < 6; ++ae) {
            C8charg_0 += calculateContribution(h50, 3, 3, ie, ae, false) + Mch / mb_muW * calculateContribution(h60, 3, 4, ie, ae, false);
            C8_chargeps_0 += Mch / mb_muW * calculateContribution(h60, 3, 4, ie, ae, true);
        }
    }

    complex_t C8H_0 = 1. / 3. * lu * lu * F8_1(yt) - lu * ld * F8_2(yt);
    
    return scalar_t(C8SMeps_0 + C8Heps_0 + C8Heps2_0 + C8charg_0 + C8_chargeps_0 + C8H_0);
}

scalar_t C8_susy::compute_NLO(const ParamSrc& src) {
    // double xt      = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1});
    double lu      = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 7);
    double ld      = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 8);
    double yt      = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    double Q_match = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);;
    double mW      = src.get_val(ParameterType::SM, "MASS", 24);
    double mH      = src.get_val(ParameterType::BSM, "MASS", 37);
    double mb_muW  = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1});

    complex_t C8charg_1 = 0.0;
    complex_t C8four_1  = 0.0;

    for (int ie = 0; ie < 2; ++ie) {
        double Mchi = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie});
        double ratio_mW_Mchi2 = std::pow(mW / Mchi, 2.0);

        for (int ae = 0; ae < 6; ++ae) {
            double MsqU = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae});
            double x = std::pow(MsqU / Mchi, 2.0);
            double log_mu = std::log(std::pow(Q_match / MsqU, 2.0));

            double A1 = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1});
            double A2 = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2});
            double B2 = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2});

            C8charg_1 += ratio_mW_Mchi2 * (
                A1 * A2 * h51(x, log_mu) +
                Mchi / mb_muW * A1 * B2 * h61(x, log_mu)
            );

            for (int ce = 0; ce < 6; ++ce) {
                double MsqU_ce = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce});
                double log_ce = std::log(std::pow(Q_match / MsqU, 2.0));
                double ratio_ce_Mchi = MsqU_ce / Mchi;

                for (int de = 0; de < 6; ++de) {
                    double x_ae = std::pow(MsqU / Mchi, 2.0);
                    double x_de = pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}) / Mchi, 2.0);

                    double Nij = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,ae, ce});
                    double Ncd = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,ce, de});
                    double A3 = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1});
                    double B3 = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, de, 2});
                    double B4 = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4,ie, de, 2});

                    C8four_1 += ratio_mW_Mchi2 * Nij * ratio_ce_Mchi * Ncd * (1.0 + log_ce) * (
                        A3 * B3 * q21(x_ae, x_de) +
                        Mchi / mb_muW * A3 * B4 * q41(x_ae, x_de)
                    );
                }
            }
        }
    }
    C8charg_1 *= -0.5 * src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 19);
    C8four_1 *= -0.5 * src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 19);
    complex_t C8H_1 = G8H(yt, lu, ld) + Delta8H(yt, lu, ld) * std::log(std::pow(Q_match / mH, 2.0)) - 1. / 6. * EH(yt, lu);

    // printf("C8H_1 (NLO) : %.8lf\n", C8H_1);
    // printf("C8charg_1 (NLO) : %.8lf\n", C8charg_1);
    // printf("C8four_1 (NLO) : %.8lf\n", C8four_1);

    return scalar_t(C8charg_1 + C8four_1 + C8H_1);
}

scalar_t C8_susy::compute_NNLO(const ParamSrc& src) {
    double lu            = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 7);
    double ld            = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 8);
    double yt            = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    double Q_match       = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);;
    double mass_top_muW  = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 6);

    double log_ratio = std::log(std::pow(Q_match / mass_top_muW, 2.0));

    return scalar_t(C8H2(yt, lu, ld, log_ratio));
}


C9_susy::C9_susy() : WilsonCoefficient("C9_SUSY", GroupMapper::str(WGroup::B, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_SI_BSM", 7},        // lu
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},     // yt
            {ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}},  // xt
            {ParameterType::WILSON, "WPARAM_SI_SM", 4},         // sw2
            {ParameterType::WILSON, "MATRIX_BSM", 13},          // B90c
            {ParameterType::WILSON, "MATRIX_BSM", 14},          // C90c
            {ParameterType::WILSON, "MATRIX_BSM", 15}           // D90c
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };

    matching_info[QCDOrder::NLO] = {
        {
            {ParameterType::WILSON, "WPARAM_SI_BSM", 19},   // kappa
            {ParameterType::WILSON, "WPARAM_SI_BSM", 7},   // lu
            // {ParameterType::WILSON, "WPARAM_SI_BSM", 8},   // ld
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", 1}, // yt
            {ParameterType::WILSON, "EW_SCALE", 1},         // Q_match
            {ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}}, // xt
            {ParameterType::SM, "MASS", 24},                // mW
            {ParameterType::BSM, "MASS", 37},               // mH
            {ParameterType::WILSON, "WPARAM_SI_SM", 4},     // sw2
            {ParameterType::SM, "GAUGE", 2},                // g2
        
            // Pour les boucles : masses
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,2}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,3}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,4}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,5}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {16,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {16,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {16,2}},
        
            
            {ParameterType::WILSON, "MATRIX_BSM", {6,0,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {6,1,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {6,0,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {6,1,1,1}},
        
            
            // Mixages
            {ParameterType::BSM, "UMIX", {0+1, 0+1}},
            {ParameterType::BSM, "UMIX", {1+1,0+1}},
            {ParameterType::BSM, "VMIX", {0+1, 0+1}},
            {ParameterType::BSM, "VMIX", {1+1,0+1}}
        },
        compute_NLO,
        get_lhaid_from_name(QCDOrder::NLO)
    };

    auto& sources = matching_info[QCDOrder::NLO].sources;

    for (int me = 0; me < 6; ++me) {
        for (int be = 0; be < 3; ++be) {
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {1, me, be}});
        }
    }

    for (int ie = 0; ie < 2; ie++) {
        for (int ae = 0; ae < 6; ae++) {
            for (int be = 1; be<3; ++be) {
                sources.insert({ParameterType::WILSON, "MATRIX_BSM", {3, ie, ae, be}});
            }
        }
    }

    for (int je = 0; je < 2; je++) {
        for (int be=0; be<3; be++) {
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {5,je, be, 1}});
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {6,je, be, 1}});
        }
    }

    for (int ce = 0; ce < 6; ce++) {
        for (int de = 0; de < 6; de++) {
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {9,ce, de}});
        }
    }
}

scalar_t C9_susy::compute_LO(const ParamSrc& src) {
    double xt   = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1});
    double lu   = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 7);
    double yt   = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    double sw2  = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);

    complex_t B90c = src.get_val(ParameterType::WILSON, "MATRIX_BSM", 13);
    complex_t C90c = src.get_val(ParameterType::WILSON, "MATRIX_BSM", 14);
    complex_t D90c = src.get_val(ParameterType::WILSON, "MATRIX_BSM", 15);

    // std::cout << "B90c " <<  B90c << std::endl;
    // std::cout << "C90c " <<  C90c << std::endl;
    // std::cout << "D90c " <<  D90c << std::endl;

    complex_t C9charg_0 = (1.0 - 4.0 * sw2) / sw2 * C90c - B90c / sw2 - D90c;
    double    C9H_0     = (1.0 - 4.0 * sw2) / sw2 * C9llH0(xt, yt, lu) - D9H0(yt, lu);

    // std::cout << "C9H_0" << C9H_0 << std::endl;
    // std::cout << "C9charg_0" << C9charg_0 << std::endl;
    return scalar_t(C9charg_0 + C9H_0);
}

scalar_t C9_susy::compute_NLO(const ParamSrc& src) {
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1});

    double lu = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 7);
    // double ld = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 8);
    double yt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    double Q_match = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);;
    double mW = src.get_val(ParameterType::SM, "MASS", 24);
    double mH = src.get_val(ParameterType::BSM, "MASS", 37);
    double sw2 = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);
    double g2 = src.get_val(ParameterType::SM, "GAUGE", 2);
    complex_t C91f = 0.0;
    complex_t B1f1 = 0.0;
    complex_t B1f2 = 0.0;
    complex_t D91f = 0.0;
    complex_t B1c1=0.;
    complex_t B1c2=0.;
    complex_t D91c= 0.;
    complex_t C91c = 0.;

    for (int ie = 0; ie < 2; ie++) {

        for (int ae = 0; ae < 6; ae++) {
            double mass24_Mch_ie_squared = pow(mW / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0);
            double ratio_MsqU_ae_Mch_ie = pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0);
            double log_mu_W_MsqU_ae = log(pow(Q_match / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}), 2.0));

            D91c += std::pow(mW / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}) *
                h31(ratio_MsqU_ae_Mch_ie, log_mu_W_MsqU_ae);

            

            for (int je = 0; je < 2; je++) {
                double ratio_Mch_je_ie = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie});
                double factor_abs = 2.0 * std::fabs(ratio_Mch_je_ie);
                double factor_f31_f30 = (f31(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) +
                                        4.0 * (f30(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) + 
                                        (f30(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 1.0001) - 
                                        f30(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 0.9999)) / 0.0002) * log_mu_W_MsqU_ae);
                double factor_f41_f40 = (f41(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) +
                                        4.0 * (f40(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) + 
                                        (f40(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 1.0001) - 
                                        f40(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 0.9999)) / 0.0002) * log_mu_W_MsqU_ae);

                C91c += src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,je, ae, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}) *
                        ((factor_abs * factor_f31_f30 * src.get_val(ParameterType::BSM, "UMIX", {je+1, 0+1}) * src.get_val(ParameterType::BSM, "UMIX", {ie+1, 0+1})) -
                        (factor_f41_f40 * src.get_val(ParameterType::BSM, "VMIX", {je+1, 0+1}) * src.get_val(ParameterType::BSM, "VMIX", {ie+1, 0+1})));

                for (int de=0; de<6; de++) {
                    for (int ke=0; ke<6; ke++) {
                        C91f+= src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,de, ke}) * pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.) *  src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,ke, ae}) * (1+log(pow(Q_match/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}),2.)))*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,je, de, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}) * (
                            2.*fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})) * f60(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.), pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.), pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.))*src.get_val(ParameterType::BSM, "UMIX", {je+1, 0+1})*src.get_val(ParameterType::BSM, "UMIX", {ie+1, 0+1})
                            -f50(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.), pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.), pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.)) *  src.get_val(ParameterType::BSM, "VMIX", {je+1, 0+1})*src.get_val(ParameterType::BSM, "VMIX", {ie+1, 0+1}));
                        
                    }
                }

                for (int be=0; be<3; be++) {
                    double ratio_Mch = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie});
                    double ratio_MsqU = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie});
                    double ratio_Msn = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie});
                    
                    B1c1 += src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,je, ae, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}) / pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2) * 
                            (0.5 * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5, ie, be, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5, je, be, 1}) * 
                            (f81(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) + 
                            4 * (f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) +
                            (f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 1.0001, pow(ratio_Msn, 2)) - 
                            f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 0.9999, pow(ratio_Msn, 2))) / 0.0002) *
                            log(pow(Q_match / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}), 2))));

                    B1c2 += src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,je, ae, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}) / pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2) * 
                            (src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6, ie, be, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6, je, be, 1}) * fabs(ratio_Mch) * 
                            (f91(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) + 
                            4 * (f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) +
                            (f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 1.0001, pow(ratio_Msn, 2)) - 
                            f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 0.9999, pow(ratio_Msn, 2))) / 0.0002) *
                            log(pow(Q_match / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}), 2))));
                }
            }
            for (int ce = 0; ce < 6; ce++) {
                for (int fe = 0; fe < 3; fe++) {
                    for (int de = 0; de < 6; de++) {
                        for (int ke = 0; ke < 6; ke++) {
                            double MsqU_ke_Mch_ie_squared = pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0);
                            double log_scale_MsqU_ke = log(pow(Q_match / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}), 2.0));
                            C91f +=  src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,de, ke}) * MsqU_ke_Mch_ie_squared *  src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,ke, ae}) * 
                                    (1.0 + log_scale_MsqU_ke) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ce, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}) * 
                                    f50(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                        pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                        pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0)) * 
                                    src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ce, fe}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, de, fe});

                            C91f +=  src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,de, ke}) * MsqU_ke_Mch_ie_squared *  src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,ke, ce}) * 
                                    (1.0 + log_scale_MsqU_ke) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, de, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}) * 
                                    f50(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                        pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                        pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0)) * 
                                    src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ce, fe}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe});
                            
                        }

                        for (int je = 0; je < 2; je++) {
                            double factor_common = mass24_Mch_ie_squared * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,ae, de}) *pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2)*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,de, ce}) *
                            (1+log(pow(Q_match/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}),2))) *  src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,je, ae, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ce, 2});

                            B1f1 += factor_common * 0.5 * f90(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                        pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                        pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                        pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, fe}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0)) *
                                    src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5, ie, fe, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5, je, fe, 1});

                            B1f2 += factor_common * fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})) * 
                                    f100(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                            pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                            pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                            pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, fe}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0)) *
                                    src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6, ie, fe, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6, je, fe, 1});
                        }
                    }
                    C91c += src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ce, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}) *
                        (f51(ratio_MsqU_ae_Mch_ie, pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0)) +
                            4.0 * (f40(ratio_MsqU_ae_Mch_ie, pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0)) + 
                            (f40(ratio_MsqU_ae_Mch_ie, pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0)* 1.0001) - 
                            f40(ratio_MsqU_ae_Mch_ie, pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0)* 0.9999)) / 0.0002
                            +(f40(ratio_MsqU_ae_Mch_ie*1.0001, pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))
                            -f40(ratio_MsqU_ae_Mch_ie*0.9999, pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0)))/0.0002)* log_mu_W_MsqU_ae) *
                        src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ce, fe}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe});
                }

                double MsqU_ce_Mch_ie_squared = pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0);
                double log_scale_MsqU_ce = 1.0 + log(pow(Q_match / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}), 2.0));

                
                for (int de = 0; de < 6; de++) {
                    // double ratio_MsqU_ae_Mch_ie = pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0);
                    double ratio_MsqU_de_Mch_ie = pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0);
                    
                    D91f += mass24_Mch_ie_squared * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,ae, ce}) * MsqU_ce_Mch_ie_squared * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,ce, de}) * log_scale_MsqU_ce * 
                        src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, de, 2}) * 
                        q51(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), ratio_MsqU_de_Mch_ie);

                }
            }
        }
    }


    double kappa = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 19);
    C91c *= -kappa / 8.0;
    D91c *= kappa;
    complex_t B91c = -(B1c1 - B1c2) * kappa * mW * mW / 2.0 / pow(g2, 2);

    C91f *= kappa / 6.0;
    D91f *= kappa;

    complex_t B91f = (B1f1 - B1f2) * 2.0 / 3.0 * kappa / pow(g2, 2);

    complex_t C9four_1 = (1. - 4. * sw2) / sw2 * C91f - B91f / sw2 - D91f;

    complex_t C9charg_1=(1.-4.*sw2)/sw2*C91c-B91c/sw2-D91c;
    complex_t C9H_1 = (1.-4.*sw2)/sw2*C9llH1(xt,yt,lu,log(pow(Q_match/mH,2.)))
    -D9H1(yt,lu,log(pow(Q_match/mH,2.)));

    return scalar_t(C9four_1+C9charg_1 + C9H_1);
    
}

C10_susy::C10_susy() : WilsonCoefficient("C10_SUSY", GroupMapper::str(WGroup::B, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_SI_BSM", 7},        // lu
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},     // yt
            {ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}},  // xt
            {ParameterType::WILSON, "WPARAM_SI_SM", 4},         // sw2
            {ParameterType::WILSON, "MATRIX_BSM", 14},          // C90c
            {ParameterType::WILSON, "MATRIX_BSM", 16}           // B100c
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };

    matching_info[QCDOrder::NLO] = {
        {
            {ParameterType::WILSON, "WPARAM_SI_BSM", 19},
            {ParameterType::WILSON, "WPARAM_SI_BSM", 7},
            // {ParameterType::WILSON, "WPARAM_SI_BSM", 8},
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},
            {ParameterType::WILSON, "EW_SCALE", 1},
            {ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}},
            {ParameterType::SM, "MASS", 24},
            {ParameterType::BSM, "GAUGE", 2},
            {ParameterType::WILSON, "WPARAM_SI_SM", 4},
            {ParameterType::BSM, "MASS", 37},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,2}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,3}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,4}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,5}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {16,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {16,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {16,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,0,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,1,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,2,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,3,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,3,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,4,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,4,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,5,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,5,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,5,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,3,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,4,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,5,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,5,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,3,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,4,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,5,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,5,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {5,0,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {5,0,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {5,0,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {5,1,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {5,1,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {5,1,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {6,0,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {6,0,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {6,0,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {6,1,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {6,1,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {6,1,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,0,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,0,3}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,0,4}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,0,5}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,1,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,1,3}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,1,4}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,1,5}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,2,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,2,3}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,2,4}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,2,5}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,3,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,3,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,3,3}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,3,4}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,3,5}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,4,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,4,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,4,3}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,4,4}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,4,5}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,5,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,5,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,5,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,5,3}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,5,4}},
            {ParameterType::WILSON, "MATRIX_BSM", {9,5,5}},
            {ParameterType::BSM, "UMIX", {1, 1}},
            {ParameterType::BSM, "UMIX", {1, 2}},
            {ParameterType::BSM, "UMIX", {2, 1}},
            {ParameterType::BSM, "UMIX", {2, 2}},
            {ParameterType::BSM, "VMIX", {1, 1}},
            {ParameterType::BSM, "VMIX", {1, 2}},
            {ParameterType::BSM, "VMIX", {2, 1}},
            {ParameterType::BSM, "VMIX", {2, 2}},
        },
        compute_NLO,
        get_lhaid_from_name(QCDOrder::NLO)
    };

}

scalar_t C10_susy::compute_LO(const ParamSrc& src) {
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1});
    double lu = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 7);
    double yt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    double sw2 = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);

    complex_t C90c   = src.get_val(ParameterType::WILSON, "MATRIX_BSM", 14);
    complex_t B100c  = src.get_val(ParameterType::WILSON, "MATRIX_BSM", 16);

    complex_t C10charg_0 = (B100c - C90c) / sw2;
    complex_t C10H_0     = -C9llH0(xt, yt, lu) / sw2;

    return C10charg_0 + C10H_0;
}

scalar_t C10_susy::compute_NLO(const ParamSrc& src) {
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1});
    double lu = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 7);
    // double ld = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 8);
    double yt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    double Q_match = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);;
    double mW = src.get_val(ParameterType::SM, "MASS", 24);
    double mH = src.get_val(ParameterType::BSM, "MASS", 37);
    double sw2 = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);
    double g2 = src.get_val(ParameterType::BSM, "GAUGE", 2);

    complex_t C91f = 0.0;
    complex_t B1f1 = 0.0;
    complex_t B1f2 = 0.0;
    complex_t B1c1=0.;
    complex_t B1c2=0.;
    complex_t C91c = 0.;
    complex_t C91f_b = 0.;
    complex_t C91f_50 = 0.;
    complex_t C91f_60 = 0.;
    complex_t C91f_beg = 0.;
    complex_t C91f_beg_1 = 0;
    complex_t C91f_beg_1a = 0;
    complex_t C91f_beg_1b = 0;
    complex_t C91f_beg_10 = 0;
    complex_t C91f_beg_11 = 0;
    complex_t C91f_a1 = 0.;
    complex_t C91f_a2 = 0.;
    // std::cout << "Qmatch : " << Q_match << std::endl;
    long long aa = 0;
    for (int ie = 0; ie < 2; ie++) {
        for (int ae = 0; ae < 6; ae++) {
            double mass24_Mch_ie_squared = pow(mW / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0);
            double ratio_MsqU_ae_Mch_ie = pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0);
            double log_mu_W_MsqU_ae = log(pow(Q_match / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}), 2.0));

            for (int je = 0; je < 2; je++) {
                double ratio_Mch_je_ie = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie});
                double factor_abs = 2.0 * std::fabs(ratio_Mch_je_ie);
                double factor_f31_f30 = (f31(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) +
                                        4.0 * (f30(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) + 
                                        (f30(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 1.0001) - 
                                        f30(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 0.9999)) / 0.0002) * log_mu_W_MsqU_ae);
                double factor_f41_f40 = (f41(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) +
                                        4.0 * (f40(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) + 
                                        (f40(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 1.0001) - 
                                        f40(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 0.9999)) / 0.0002) * log_mu_W_MsqU_ae);

                C91c += src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,je, ae, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}) *
                        ((factor_abs * factor_f31_f30 * src.get_val(ParameterType::BSM, "UMIX", {je+1,0+1}) * src.get_val(ParameterType::BSM, "UMIX", {ie+1, 0+1})) -
                        (factor_f41_f40 * src.get_val(ParameterType::BSM, "VMIX", {je+1,0+1}) * src.get_val(ParameterType::BSM, "VMIX", {ie+1, 0+1})));


                for (int de=0; de<6; de++) {
                    for (int ke=0; ke<6; ke++) {
                        C91f+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,de, ke}) * pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,ke, ae}) * (1.+log(pow(Q_match/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}),2.)))*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,je, de, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}) * (
                            2.*fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})) * f60(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.), pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.), pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.))*src.get_val(ParameterType::BSM, "UMIX", {je+1,0+1})*src.get_val(ParameterType::BSM, "UMIX", {ie+1,0+1})
                            -f50(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.), pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.), pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.)) *  src.get_val(ParameterType::BSM, "VMIX", {je+1,0+1})*src.get_val(ParameterType::BSM, "VMIX", {ie+1,0+1}));
                        
                            C91f_b +=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,de, ke}) * pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,ke, ae}) * (1+log(pow(Q_match/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}),2.)))*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,je, de, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}) * (
                            2.*fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})) * f60(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.), pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.), pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.))*src.get_val(ParameterType::BSM, "UMIX", {je+1,0+1})*src.get_val(ParameterType::BSM, "UMIX", {ie+1,0+1})
                            -f50(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.), pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.), pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.)) *  src.get_val(ParameterType::BSM, "VMIX", {je+1,0+1})*src.get_val(ParameterType::BSM, "VMIX", {ie+1,0+1}));
                            
                            C91f_60 +=  2.*fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})) * f60(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.), pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.), pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.))*src.get_val(ParameterType::BSM, "UMIX", {je+1,0+1})*src.get_val(ParameterType::BSM, "UMIX", {ie+1,0+1});
                            C91f_50 += f50(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.), pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.), pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.)) *  src.get_val(ParameterType::BSM, "VMIX", {je+1,0+1})*src.get_val(ParameterType::BSM, "VMIX", {ie+1,0+1});
                            C91f_beg += src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,de, ke}) * pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,ke, ae}) * (1+log(pow(Q_match/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}),2.)))*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,je, de, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2});
                            
                            C91f_beg_1 +=  (1.+log(pow(Q_match/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}),2.)))*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,je, de, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2});
                            C91f_beg_1a +=  (1.+log(pow(Q_match/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}),2.)))*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,je, de, 1});
                            C91f_beg_1b +=  (1.+log(pow(Q_match/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}),2.)))*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2});
                            // C91f_beg_12 +=  (1+log(pow(Q_match/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}),2.)))*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,je, de, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2});
                            C91f_beg_11 +=  src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,je, de, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2});
                            C91f_beg_10 +=  (1.+log(pow(Q_match/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}),2.)));
                        aa++;        
                        }
                }
                
                for (int be=0; be<3; be++) {
                    double ratio_Mch = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie});
                    double ratio_MsqU = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie});
                    double ratio_Msn = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie});
                    
                    B1c1 += src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,je, ae, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}) / pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2) * 
                            (0.5 * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5, ie, be, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5, je, be, 1}) * 
                            (f81(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) + 
                            4 * (f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) +
                            (f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 1.0001, pow(ratio_Msn, 2)) - 
                            f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 0.9999, pow(ratio_Msn, 2))) / 0.0002) *
                            log(pow(Q_match / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}), 2))));

                    B1c2 += src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,je, ae, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}) / pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2) * 
                            (src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6, ie, be, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6, je, be, 1}) * fabs(ratio_Mch) * 
                            (f91(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) + 
                            4 * (f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) +
                            (f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 1.0001, pow(ratio_Msn, 2)) - 
                            f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 0.9999, pow(ratio_Msn, 2))) / 0.0002) *
                            log(pow(Q_match / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}), 2))));
                }
            }

            for (int ce = 0; ce < 6; ce++) {
                for (int fe = 0; fe < 3; fe++) {
                    for (int de = 0; de < 6; de++) {
                        for (int ke = 0; ke < 6; ke++) {
                            double MsqU_ke_Mch_ie_squared = pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0);
                            double log_scale_MsqU_ke = log(pow(Q_match / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}), 2.0));
                            C91f += src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,de, ke}) * MsqU_ke_Mch_ie_squared * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,ke, ae}) * 
                                    (1.0 + log_scale_MsqU_ke) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ce, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}) * 
                                    f50(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                        pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                        pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0)) * 
                                    src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ce, fe}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, de, fe});

                            C91f += src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,de, ke}) * MsqU_ke_Mch_ie_squared * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,ke, ce}) * 
                                    (1.0 + log_scale_MsqU_ke) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, de, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}) * 
                                    f50(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                        pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                        pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0)) * 
                                    src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ce, fe}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe});
                            
                            C91f_a1 += src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,de, ke}) * MsqU_ke_Mch_ie_squared * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,ke, ae}) * 
                                    (1.0 + log_scale_MsqU_ke) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ce, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}) * 
                                    f50(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                        pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                        pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0)) * 
                                    src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ce, fe}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, de, fe});

                            C91f_a2 += src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,de, ke}) * MsqU_ke_Mch_ie_squared * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,ke, ce}) * 
                                    (1.0 + log_scale_MsqU_ke) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, de, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}) * 
                                    f50(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                        pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                        pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0)) * 
                                    src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ce, fe}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe});

                        }

                        for (int je = 0; je < 2; je++) {
                            double factor_common = mass24_Mch_ie_squared * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,ae, de}) *pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2)*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,de, ce}) *
                            (1+log(pow(Q_match/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}),2))) *  src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,je, ae, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ce, 2});

                            B1f1 += factor_common * 0.5 * f90(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                        pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                        pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                        pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, fe}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0)) *
                                    src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5, ie, fe, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5, je, fe, 1});

                            B1f2 += factor_common * fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})) * 
                                    f100(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                            pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                            pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0), 
                                            pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, fe}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0)) *
                                    src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6, ie, fe, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6, je, fe, 1});

                            
                        }
                    }

                    C91c += src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ce, 1}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}) *
                        (f51(ratio_MsqU_ae_Mch_ie, pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0)) +
                            4.0 * (f40(ratio_MsqU_ae_Mch_ie, pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0)) + 
                            (f40(ratio_MsqU_ae_Mch_ie, pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0)* 1.0001) - 
                            f40(ratio_MsqU_ae_Mch_ie, pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0)* 0.9999)) / 0.0002
                            +(f40(ratio_MsqU_ae_Mch_ie*1.0001, pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))
                            -f40(ratio_MsqU_ae_Mch_ie*0.9999, pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}) / src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2.0)))/0.0002)* log_mu_W_MsqU_ae) *
                        src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ce, fe}) * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe});
                }

            }
        }
    }

    // for(int ie=0;ie<2;ie++) for(int ae=0;ae<6;ae++) printf("X_UL[%d][%d][1] : %lf\n",ie, ae, src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1}).real());
    // printf("aa : %d\n", aa);
    // printf("C91f before : %lf\n", C91f_b.real());
    // printf("C91f 50 : %lf\n", C91f_50.real());
    // printf("C91f 60 : %lf\n", C91f_60.real());
    // printf("C91f beg : %lf\n", C91f_beg.real());
    // printf("C91f beg1 : %lf\n", C91f_beg_1.real());
    // printf("C91f beg1a : %lf\n", C91f_beg_1a.real());
    // printf("C91f beg1b : %lf\n", C91f_beg_1b.real());
    // printf("C91f beg10 : %lf\n", C91f_beg_10.real());
    // printf("C91f beg11 : %lf\n", C91f_beg_11.real());
    // printf("C91f beg10 img : %lf\n", C91f_beg_10.imag());
    // printf("C91f beg11 img : %lf\n", C91f_beg_11.imag());
    // printf("C91f a1 : %lf\n", C91f_a1.real());
    // printf("C91f a1 : %lf\n", C91f_a2.real());
    // for (int ie =0; ie<2; ++ie) {
    //     printf("U_MIX[%d][1] = %lf\n", ie+1, src.get_val(ParameterType::BSM, "UMIX", {ie+1,0+1}).real());
    // }
    // for (int ie =0; ie<2; ++ie) {
    //     printf("V_MIX[%d][1] = %lf\n", ie+1, src.get_val(ParameterType::BSM, "VMIX", {ie+1,0+1}).real());
    // }
    double kappa = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 19);
    C91c *= -kappa / 8.0;

    complex_t B101c = (B1c1 + B1c2) * kappa * mW * mW / 2.0 / pow(g2, 2);

    C91f *= kappa / 6.0;

    complex_t B101f = -(B1f1 + B1f2) * 2.0 / 3.0 * kappa / pow(g2, 2);

    // printf("B101f : %f\n", B101f);
    // printf("C91f : %f\n", C91f);
    // printf("sw2 : %f\n", sw2);

    complex_t C10four_1 = (B101f - C91f) / sw2;	

    complex_t C10charg_1=(B101c-C91c)/sw2;
    double C10H_1 = -C9llH1(xt,yt,lu,log(pow(Q_match/mH,2.)))/sw2;

    // printf("C10H_1 : %f\n", C10H_1);
    // printf("C10charg_1 : %f\n", C10charg_1);
    // printf("C10four_1 : %f\n", C10four_1);

    return C10four_1 + C10charg_1 + C10H_1;
}


CP7_susy::CP7_susy() : WilsonCoefficient("CP7_SUSY", GroupMapper::str(WGroup::BPrime, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_SI_BSM", 19},
            {ParameterType::WILSON, "WPARAM_SI_BSM", 8},
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},
            {ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}},
            {ParameterType::SM, "MASS", 24},
            {ParameterType::SM, "MASS", 3},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13, 0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13, 1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 2}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 3}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 4}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14, 5}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 0, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 0, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 1, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 1, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 2, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 2, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 3, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 3, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 4, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 4, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 5, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 0, 5, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 0, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 0, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 1, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 1, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 2, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 2, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 3, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 3, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 4, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 4, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 5, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4, 1, 5, 2}},
            
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 0, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 0, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 1, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 1, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 2, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 2, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 3, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 3, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 4, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 4, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 5, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 0, 5, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 0, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 0, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 1, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 1, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 2, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 2, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 3, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 3, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 4, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 4, 2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 5, 1}},
            {ParameterType::WILSON, "MATRIX_BSM", {3, 1, 5, 2}}
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };
}

scalar_t CP7_susy::compute_LO(const ParamSrc& src) {
    double ld = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 8);
    double yt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    double mass_top_muW = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 6);
    double mass_b_muW = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1});
    double mW = src.get_val(ParameterType::SM, "MASS", 24);
    double C7pH = src.get_val(ParameterType::SM, "MASS", 3) * mass_b_muW / pow(mass_top_muW, 2) * (1.0 / 3.0) * ld * ld * F7_1(yt);

    double C7pcharg = 0.0;
    for (int ie = 0; ie < 2; ++ie) {
        for (int ae = 0; ae < 6; ++ae) {
            double Mch = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie});
            double MsqU = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae});
            double ratio = pow(MsqU / Mch, 2);
            C7pcharg += pow(mW / Mch, 2) *
                (src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 1}) *
                 src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2}) * h10(ratio)
                + Mch / mass_b_muW *
                  src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 1}) *
                  src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, ie, ae, 2}) * h20(ratio));
        }
    }

    double kappa = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 19);
    C7pcharg *= -0.5 * kappa;

    return C7pH + C7pcharg;
}

CP8_susy::CP8_susy() : WilsonCoefficient("CP8_SUSY", GroupMapper::str(WGroup::BPrime, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_SI_BSM", 19},
            {ParameterType::WILSON, "WPARAM_SI_BSM", 8},
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},
            {ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}},
            {ParameterType::SM, "MASS", 24},
            {ParameterType::SM, "MASS", 3},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,2}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,3}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,4}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,5}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,3,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,4,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,5,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,5,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,0,5,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,3,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,4,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,5,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,5,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {3,1,5,2}},
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };
}

scalar_t CP8_susy::compute_LO(const ParamSrc& src) {
    double ld = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 8);
    double yt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    double mt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 6);
    double mb = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1});
    double mW = src.get_val(ParameterType::SM, "MASS", 24);
    double m3 = src.get_val(ParameterType::SM, "MASS", 3);

    double C8pH = m3 * mb / (mt * mt) * (1.0 / 3.0) * ld * ld * F8_1(yt);

    double C8pcharg = 0.0;
    for (int ie = 0; ie < 2; ++ie) {
        double Mch = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie});
        for (int ae = 0; ae < 6; ++ae) {
            double MsqU = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae});
            double x = pow(MsqU / Mch, 2);

            C8pcharg += pow(mW / Mch, 2) * (
                src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 1}) *
                src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2}) * h50(x) +
                Mch / mb *
                src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 1}) *
                src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}) * h60(x));
        }
    }

    double kappa = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 19);
    C8pcharg *= -0.5 * kappa;

    return C8pH + C8pcharg;
}

CP9_susy::CP9_susy() : WilsonCoefficient("CP9_SUSY", GroupMapper::str(WGroup::BPrime, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_SI_BSM", 19},   // kappa
            // {ParameterType::WILSON, "WPARAM_SI_BSM", 7},   // lu
            {ParameterType::WILSON, "WPARAM_SI_BSM", 8},   // ld
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},   // yt
            // {ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}},  // xt
            {ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}},  // mass_b_muW
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},  // mass_top_muW
            {ParameterType::WILSON, "WPARAM_SI_SM", 3},     // used by C10pH
            {ParameterType::WILSON, "WPARAM_SI_SM", 4},     // sw2
            {ParameterType::SM, "MASS", 3},                 // mass used by C9pH and C10pH
            {ParameterType::SM, "MASS", 24},                // mW
            {ParameterType::SM, "GAUGE", 2},                // g2
            {ParameterType::BSM, "MASS", 37},               // mH
            {ParameterType::BSM, "MINPAR", 3},                // tanb
            // {ParameterType::WILSON, "EW_SCALE", 1},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,2}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,3}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,4}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,5}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {16,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {16,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {16,2}},
        
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,3,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,4,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,5,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,5,2}},
        
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,3,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,4,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,5,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,5,2}},
        
            {ParameterType::WILSON, "MATRIX_BSM", {5,0,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {5,0,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {5,0,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {5,1,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {5,1,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {5,1,2,1}},
        
            {ParameterType::WILSON, "MATRIX_BSM", {6,0,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {6,0,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {6,0,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {6,1,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {6,1,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {6,1,2,1}},
        
            {ParameterType::WILSON, "MATRIX_BSM", {2,0,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,1,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,2,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,3,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,3,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,4,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,4,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,5,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,5,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,5,2}},
        
            {ParameterType::BSM, "UMIX", {0+1, 0+1}},
            {ParameterType::BSM, "UMIX", {1+1,0+1}},
            {ParameterType::BSM, "VMIX", {0+1, 0+1}},
            {ParameterType::BSM, "VMIX", {1+1,0+1}}
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };
}

scalar_t CP9_susy::compute_LO(const ParamSrc& src) {
    // double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1});

    // double lu = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 7);
    double ld = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 8);
    double yt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    // double Q_match = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);;
    double mW = src.get_val(ParameterType::SM, "MASS", 24);
    double mH = src.get_val(ParameterType::BSM, "MASS", 37);
    double sw2 = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);
    double mass_b_muW = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1});
    double mass_top_muW = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 6);
    double g2 = src.get_val(ParameterType::SM, "GAUGE", 2);
    double tanb = src.get_val(ParameterType::BSM, "MINPAR", 3);

    double B10pc=0.;
    double C9pc=0.;
    double B9pc=0.;
    double D9pc=0.; 

    for(int ie=0;ie<2;ie++) {
        for(int ae=0;ae<6;ae++) {
            for(int je=0;je<2;je++) {
                C9pc+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4,je, ae, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2})*(2.*fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}))*f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))*src.get_val(ParameterType::BSM, "VMIX", {je+1, 0+1})*src.get_val(ParameterType::BSM, "VMIX", {ie+1, 0+1}) -f40(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))*src.get_val(ParameterType::BSM, "UMIX", {je+1, 0+1})*src.get_val(ParameterType::BSM, "UMIX", {ie+1, 0+1}));
                for(int be=0;be<6;be++) {
                    if (be<3) {
                    B10pc+=-src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4,je, ae, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*(0.5*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6,ie, be, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6,je, be, 1})*f50(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)) +src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5,ie, be, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5,je, be, 1})*fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}))*f60(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)));
                    B9pc+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4,je, ae, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*(0.5*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6,ie, be, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6,je, be, 1})*f50(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)) -src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5,ie, be, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5,je, be, 1})*fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}))*f60(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)));
                    }
                }

            }
            for(int be=0;be<6;be++) {
                for(int ce=0;ce<3;ce++) {
                    C9pc+=-src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4,ie, be, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2})*f40(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2,be, ce})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2,ae, ce});
                    }
            }
            D9pc+=pow(mW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2})*h30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.));
        }
    } 		
    
    double kappa = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 19);

    B9pc*=kappa*(mW)*mW/2./g2/g2;	

    C9pc*=-kappa/8.;

    
    double C10pH = -mass_b_muW*(src.get_val(ParameterType::SM, "MASS", 3))*(tanb*tanb/8./mW/mW
    +pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 3)*tanb*tanb/4./mW/mH,2.))*f20(yt)/sw2;
    
    double C9pH =(4.*sw2-1.)*C10pH - src.get_val(ParameterType::SM, "MASS", 3)*mass_b_muW/mass_top_muW/mass_top_muW*D9H0(yt,ld);
    
    
    D9pc*=kappa;

    double C9pcharg=(1.-4.*sw2)/sw2*C9pc-B9pc/sw2-D9pc;

    return C9pH+C9pcharg;
}

CP10_susy::CP10_susy() : WilsonCoefficient("CP10_SUSY", GroupMapper::str(WGroup::BPrime, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_SI_BSM", 19},
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},
            {ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}},
            {ParameterType::WILSON, "WPARAM_SI_SM", 3},
            {ParameterType::WILSON, "WPARAM_SI_SM", 4},
            {ParameterType::SM, "MASS", 24},
            {ParameterType::SM, "MASS", 3},
            {ParameterType::SM, "GAUGE", 2},
            {ParameterType::BSM, "MASS", 37},
            {ParameterType::BSM, "MINPAR", 3},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,2}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,3}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,4}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,5}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {16,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {16,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {16,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,0,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,1,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,2,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,3,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,3,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,4,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,4,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,5,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,5,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,5,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,3,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,4,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,5,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,0,5,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,3,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,3,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,4,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,4,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,5,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {4,1,5,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {5,0,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {5,0,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {5,0,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {5,1,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {5,1,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {5,1,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {6,0,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {6,0,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {6,0,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {6,1,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {6,1,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {6,1,2,1}},
            {ParameterType::BSM, "UMIX", {0+1, 0+1}},
            {ParameterType::BSM, "UMIX", {1+1, 0+1}},
            {ParameterType::BSM, "VMIX", {0+1, 0+1}},
            {ParameterType::BSM, "VMIX", {1+1, 0+1}}
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };
}

scalar_t CP10_susy::compute_LO(const ParamSrc& src) {
    double yt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    double mW = src.get_val(ParameterType::SM, "MASS", 24);
    double mH = src.get_val(ParameterType::BSM, "MASS", 37);
    double sw2 = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);
    double mass_b_muW = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1});
    double g2 = src.get_val(ParameterType::SM, "GAUGE", 2);
    double tanb = src.get_val(ParameterType::BSM, "MINPAR", 3);

    double B10pc = 0.0;
    double C9pc = 0.0;
    // printf("g2 : %.10lf\n", g2);
    // printf("mass_b_muW : %.10lf\n", mass_b_muW);
    // printf("yt : %.10lf\n", yt);
    for (int ie = 0; ie < 2; ie++) {
        for (int ae = 0; ae < 6; ae++) {
            for (int je = 0; je < 2; je++) {
                double ratio_je_ie = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}) /
                                     src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie});
                double ratio_sq = pow(ratio_je_ie, 2);
                double xae = pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}) /
                                 src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2);

                C9pc += src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4, je, ae, 1}) *
                        src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2}) *
                        (2.0 * fabs(ratio_je_ie) * f30(ratio_sq, xae) *
                         src.get_val(ParameterType::BSM, "VMIX", {je+1, 0+1}) *
                         src.get_val(ParameterType::BSM, "VMIX", {ie+1, 0+1}) -
                         f40(ratio_sq, xae) *
                         src.get_val(ParameterType::BSM, "UMIX", {je+1, 0+1}) *
                         src.get_val(ParameterType::BSM, "UMIX", {ie+1, 0+1}));

                for (int be = 0; be < 6; be++) {
                    if (be < 3) {
                        double xbe = pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}) /
                                         src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2);

                        B10pc += -src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4, je, ae, 1}) *
                                 src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2}) /
                                 pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2) *
                                 (0.5 * src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6, ie, be, 1}) *
                                  src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6, je, be, 1}) *
                                  f50(ratio_sq, xae, xbe) +
                                  src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5, ie, be, 1}) *
                                  src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5, je, be, 1}) *
                                  fabs(ratio_je_ie) *
                                  f60(ratio_sq, xae, xbe));
                    }
                }
            }
            for (int be = 0; be < 6; be++) {
                for (int ce = 0; ce < 3; ce++) {
                    double xae = pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}) /
                                     src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2);
                    double xbe = pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, be}) /
                                     src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}), 2);

                    C9pc += -src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4, ie, be, 1}) *
                             src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2}) *
                             f40(xae, xbe) *
                             src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, be, ce}) *
                             src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, ae, ce});
                }
            }
        }
    }

    double kappa = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 19);
    B10pc *= kappa * mW * mW / (2.0 * g2 * g2);
    C9pc *= -kappa / 8.0;

    // printf("B10pc : %.10lf\n", B10pc);
	// printf("C9pc : %.10lf\n", C9pc);
    // printf("kappa : %.10lf\n", kappa);
    double C10pH = -mass_b_muW * src.get_val(ParameterType::SM, "MASS", 3) *
                   (tanb * tanb / (8.0 * mW * mW) +
                    pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 3) * tanb * tanb /
                        (4.0 * mW * mH), 2)) * f20(yt) / sw2;

    double C10pcharg = (B10pc - C9pc) / sw2;

    // printf("C10pcharg : %.10lf\n", C10pcharg);
	// printf("C10pH : %.10lf\n", C10pH);
    return C10pH + C10pcharg;
}

CPQ1_susy::CPQ1_susy(WCoef coef) : WilsonCoefficient(WCoefMapper::str(coef) + "_SUSY", GroupMapper::str(WGroup::BPrime, ScaleType::MATCHING)) {
    const int lepton_mass_slot = WCoefMapper::lepton_mass_slot_from_index(WCoefMapper::lepton_index_from_cpq1(coef));
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_SI_BSM", 19},
            {ParameterType::WILSON, "WPARAM_SI_BSM", 1},
            {ParameterType::WILSON, "WPARAM_SI_BSM", 11},
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,0}},
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,1}},
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,2}},
            {ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}},
            {ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot},
            {ParameterType::WILSON, "WPARAM_SI_SM", 4},
            {ParameterType::SM, "MASS", 24},
            {ParameterType::SM, "MASS", 3},
            {ParameterType::SM, "GAUGE", 2},
            {ParameterType::BSM, "MASS", 37},
            {ParameterType::BSM, "HMIX", 1},
            {ParameterType::BSM, "MINPAR", 3},
            {ParameterType::WILSON, "EPSILON_SUSY", 5},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,2}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,3}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,4}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,5}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {16,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {16,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {16,2}},
            {ParameterType::BSM, "UMIX", {1+1, 0+1}},
            {ParameterType::BSM, "UMIX", {1+1, 1+1}},
            {ParameterType::BSM, "UMIX", {0+1, 0+1}},
            {ParameterType::BSM, "UMIX", {0+1, 1+1}},
            {ParameterType::BSM, "VMIX", {0+1, 0+1}},
            {ParameterType::BSM, "VMIX", {1+1, 1+1}},
            {ParameterType::BSM, "VMIX", {0+1, 1+1}},
            {ParameterType::BSM, "VMIX", {1+1, 0+1}}
        },
        [lepton_mass_slot](const ParamSrc& src) { return compute_LO(src, lepton_mass_slot); },
        get_lhaid_from_name(QCDOrder::LO)
    };

    auto& sources = matching_info[QCDOrder::LO].sources;

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 3; ++j) {
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {1, i, j}});
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {2, i, j}});
            // sources.insert({ParameterType::WILSON, "MATRIX_BSM", {9, i, j}});
        }
    }

    for (int ae = 0; ae < 6; ++ae) {
        for (int ie = 0; ie < 2; ++ie) {
            for (int be = 0; be < 3; ++be) {
                for (int ne = 0; ne < 3; ++ne) {
                    sources.insert({ParameterType::WILSON, "MATRIX_BSM", {12, ae, ie, be, ne}});
                }
            }
        }
    }

    for (int i = 0; i < 3; ++i) {
        sources.insert({ParameterType::WILSON, "MATRIX_BSM", {6, 0, i, 1}});
        sources.insert({ParameterType::WILSON, "MATRIX_BSM", {6, 1, i, 1}});
        sources.insert({ParameterType::WILSON, "MATRIX_BSM", {5, 0, i, 1}});
        sources.insert({ParameterType::WILSON, "MATRIX_BSM", {5, 1, i, 1}});
    }

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 6; ++j) {
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {3, i, j, 2}});
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {4, i, j, 1}});
        }
    }

}

scalar_t CPQ1_susy::compute_LO(const ParamSrc& src) {
    return compute_LO(src, WCoefMapper::lepton_mass_slot_from_index(1));
}

scalar_t CPQ1_susy::compute_LO(const ParamSrc& src, int lepton_mass_slot) {
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1});

    double z = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 1);
    double aY = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 11);
    double mW = src.get_val(ParameterType::SM, "MASS", 24);
    double mH = src.get_val(ParameterType::BSM, "MASS", 37);
    double sw2 = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);
    double g2 = src.get_val(ParameterType::SM, "GAUGE", 2);
    double tanb = src.get_val(ParameterType::BSM, "MINPAR", 3);

    double muQ = src.get_val(ParameterType::BSM, "HMIX", 1);
    double BQ1pc1=0.;
    double BQ1pc2=0.;

    double Dp{0},Dm{0};
    double a0a{0}, a0b{0}, a0c{0}, a0Q1{0}, a0Q2{0}, a1{0}, NQ1pc{0}, NQ2pc{0};

    for(int ie=0;ie<2;ie++) {
        for(int ae=0;ae<6;ae++) {

            for(int je=0;je<2;je++) {
                for(int be=0;be<6;be++) {
                    if (be<3) {
                        BQ1pc1+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4,je, ae, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})
                                /src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5,ie, be, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6,je, be, 1})
                                *f50(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})
                                /src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})
                                /src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))); 	
                        BQ1pc2+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4,je, ae, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2})
                                /src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6,ie, be, 1})
                                *src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5,je, be, 1})*fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}))
                                *f60(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})
                                /src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)));
                    }
                    for(int me=0;me<3;me++) {
                        for(int ne=0;ne<3;ne++) {
                            Dp=0.;
                            Dm=0.;
                            for(int fe=0;fe<3;fe++) { 
                                Dp+=src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe})/sqrt(2.)/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*muQ*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2,ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1,be, fe})+src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1,ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2,be, fe}));
                                Dm+=src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe})/sqrt(2.)/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*muQ*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2,ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1,be, fe})-src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1,ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2,be, fe}));
                            }
                            a0a=-(fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}))*f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.))*src.get_val(ParameterType::BSM, "UMIX", {ie+1, 1+1})*src.get_val(ParameterType::BSM, "VMIX", {je+1, 0+1}))*kron(ae,be);
                            
                            a0b=-(f40(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.))*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*src.get_val(ParameterType::BSM, "VMIX", {ie+1, 0+1}))*kron(ae,be);
                            a0c=1./mW*f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))*kron(ie,je);
                            a0Q1=a0a+a0b+Dp*a0c;
                            a0Q2=-a0a+a0b+Dm*a0c;
                            a1=src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/sqrt(2.)/mW*f80(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))*kron(ie,je)*kron(ae,be);
                            
                            NQ1pc+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, me, ne})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1,be, me})*src.get_val(ParameterType::BSM, "UMIX", {je+1, 2})*(a0Q1+a1*tanb);
                            NQ2pc+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, me, ne})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1,be, me})*src.get_val(ParameterType::BSM, "UMIX", {je+1, 2})*(a0Q2+a1*tanb);
                        }
                    }
                }

            }
        }
    }		
    
    /* Wilson coefficients CQ1 and CQ2 prime */ 
    double NQ1pH=-src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot)*(tanb)*tanb/4./mW/mW*xt*f30(xt,z);
    
    double BQ1pH=src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot)*(tanb)*tanb/4./mW/mW*f70(xt,z);
    
    complex_t CQ1pH=(NQ1pH+BQ1pH)*src.get_val(ParameterType::SM, "MASS", 3)/sw2;
    
    double kappa = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 19);
    double BQ1pc=(BQ1pc1+BQ1pc2)*kappa*(mW)*mW/2./g2/g2/sw2;

    NQ1pc*=src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot)*(tanb)*tanb/mW/(mH*mH-mW*mW)*aY*(src.get_val(ParameterType::SM, "MASS", 3))/sw2;


    complex_t CQ1pcharg=NQ1pc+BQ1pc;

    // this->set_WilsonCoeffMatching("LO", (CQ1pcharg+CQ1pH)/sus_param->epsfac);
    // return (CQ1pcharg+CQ1pH)/sus_param->epsfac;
    double epsfac = src.get_val(ParameterType::WILSON, "EPSILON_SUSY", 5);

    return (CQ1pcharg+CQ1pH)/epsfac;
}


CPQ2_susy::CPQ2_susy(WCoef coef) : WilsonCoefficient(WCoefMapper::str(coef) + "_SUSY", GroupMapper::str(WGroup::BPrime, ScaleType::MATCHING)) {
    const int lepton_mass_slot = WCoefMapper::lepton_mass_slot_from_index(WCoefMapper::lepton_index_from_cpq2(coef));
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_SI_BSM", 1},
            {ParameterType::WILSON, "WPARAM_SI_BSM", 19},
            {ParameterType::WILSON, "WPARAM_SI_BSM", 11},
            {ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot},
            {ParameterType::WILSON, "WPARAM_SI_SM", 4},
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,0}},
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,1}},
            {ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,2}},
            {ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}},
            {ParameterType::SM, "MASS", 3},
            {ParameterType::SM, "MASS", 24},
            {ParameterType::SM, "GAUGE", 2},
            {ParameterType::BSM, "MASS", 37},
            {ParameterType::BSM, "HMIX", 1},
            {ParameterType::BSM, "MINPAR", 3},
            {ParameterType::WILSON, "EPSILON_SUSY", 5},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,2}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,3}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,4}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {14,5}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {16,0}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {16,1}},
            {ParameterType::WILSON, "WPARAM_SI_BSM", {16,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,0,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,1,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,2,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {1,2,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,0,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,0,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,0,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,1,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,1,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,1,2}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,2,0}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,2,1}},
            {ParameterType::WILSON, "MATRIX_BSM", {2,2,2}},
            {ParameterType::BSM, "UMIX", {0+1, 1+1}},
            {ParameterType::BSM, "UMIX", {1+1, 1+1}},
            {ParameterType::BSM, "VMIX", {0+1, 0+1}},
            {ParameterType::BSM, "VMIX", {1+1, 0+1}}
        },
        [lepton_mass_slot](const ParamSrc& src) { return compute_LO(src, lepton_mass_slot); },
        get_lhaid_from_name(QCDOrder::LO)
    };

    auto& sources = matching_info[QCDOrder::LO].sources;

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 3; ++j) {
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {1, i, j}});
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {2, i, j}});
            // sources.insert({ParameterType::WILSON, "MATRIX_BSM", {9, i, j}});
        }
    }

    for (int ae = 0; ae < 6; ++ae) {
        for (int ie = 0; ie < 2; ++ie) {
            for (int be = 0; be < 3; ++be) {
                for (int ne = 0; ne < 3; ++ne) {
                    sources.insert({ParameterType::WILSON, "MATRIX_BSM", {12, ae, ie, be, ne}});
                }
            }
        }
    }

    for (int i = 0; i < 3; ++i) {
        sources.insert({ParameterType::WILSON, "MATRIX_BSM", {6, 0, i, 1}});
        sources.insert({ParameterType::WILSON, "MATRIX_BSM", {6, 1, i, 1}});
        sources.insert({ParameterType::WILSON, "MATRIX_BSM", {5, 0, i, 1}});
        sources.insert({ParameterType::WILSON, "MATRIX_BSM", {5, 1, i, 1}});
    }

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 6; ++j) {
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {3, i, j, 2}});
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {4, i, j, 1}});
        }
    }
}

scalar_t CPQ2_susy::compute_LO(const ParamSrc& src) {
    return compute_LO(src, WCoefMapper::lepton_mass_slot_from_index(1));
}

scalar_t CPQ2_susy::compute_LO(const ParamSrc& src, int lepton_mass_slot) {
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1});

    double z = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 1);
    double aY = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 11);
    double mW = src.get_val(ParameterType::SM, "MASS", 24);
    double mH = src.get_val(ParameterType::BSM, "MASS", 37);
    double sw2 = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);
    double g2 = src.get_val(ParameterType::SM, "GAUGE", 2);
    double tanb = src.get_val(ParameterType::BSM, "MINPAR", 3);

    double muQ = src.get_val(ParameterType::BSM, "HMIX", 1);

    double kappa = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 19);
    double epsfac = src.get_val(ParameterType::WILSON, "EPSILON_SUSY", 5);

    double BQ1pc1=0.;
    double BQ1pc2=0.;

    double Dp{0},Dm{0};
    double a0a{0}, a0b{0}, a0c{0}, a0Q1{0}, a0Q2{0}, a1{0}, NQ1pc{0}, NQ2pc{0};

    for(int ie=0;ie<2;ie++) {
        for(int ae=0;ae<6;ae++) {

            for(int je=0;je<2;je++) {

                for(int be=0;be<6;be++) {
                    if (be<3) {

                    BQ1pc1+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4,je, ae, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5,ie, be, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6,je, be, 1})*f50(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))); 	
                    BQ1pc2+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4,je, ae, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6,ie, be, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5,je, be, 1})*fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}))*f60(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)));
                    }
                    for(int me=0;me<3;me++) {
                        for(int ne=0;ne<3;ne++) {
                            Dp=0.;
                            Dm=0.;
                            for(int fe=0;fe<3;fe++) { 
                                Dp+=src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe})/sqrt(2.)/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*muQ*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, be, fe})+src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, be, fe}));
                                Dm+=src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe})/sqrt(2.)/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*muQ*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, be, fe})-src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, be, fe}));
                            }
                            a0a=-(fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}))*f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.))*src.get_val(ParameterType::BSM, "UMIX", {ie+1, 1+1})*src.get_val(ParameterType::BSM, "VMIX", {je+1, 0+1}))*kron(ae,be);

                            a0b=-(f40(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.))*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*src.get_val(ParameterType::BSM, "VMIX", {ie+1, 0+1}))*kron(ae,be);
                            a0c=1./mW*f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))*kron(ie,je);
                            a0Q1=a0a+a0b+Dp*a0c;
                            a0Q2=-a0a+a0b+Dm*a0c;
                            a1=src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/sqrt(2.)/mW*f80(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))*kron(ie,je)*kron(ae,be);
                            
                            NQ1pc+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, me, ne})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, be, me})*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*(a0Q1+a1*tanb);
                            NQ2pc+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, me, ne})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, be, me})*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*(a0Q2+a1*tanb);
                        }
                    }
                }

            }

        }
    } 		

    /* Wilson coefficients CQ1 and CQ2 prime */ 
    double NQ1pH=-src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot)*(tanb)*tanb/4./mW/mW*xt*f30(xt,z);
    
    double BQ1pH=src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot)*(tanb)*tanb/4./mW/mW*f70(xt,z);
    
    complex_t CQ1pH=(NQ1pH+BQ1pH)*src.get_val(ParameterType::SM, "MASS", 3)/sw2;
    
    complex_t CQ2pH=CQ1pH;

    
    double BQ2pc=(BQ1pc1-BQ1pc2)*kappa*(mW)*mW/2./g2/g2/sw2;

    NQ2pc*=src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot)*(tanb)*tanb/mW/(mH*mH-mW*mW)*aY*(src.get_val(ParameterType::SM, "MASS", 3))/sw2;


    complex_t CQ2pcharg=NQ2pc+BQ2pc;

    return (CQ2pH+CQ2pcharg)/epsfac;
}

CQ1_susy::CQ1_susy(WCoef coef) : WilsonCoefficient(WCoefMapper::str(coef) + "_SUSY", GroupMapper::str(WGroup::BScalar, ScaleType::MATCHING)){
    const int lepton_mass_slot = WCoefMapper::lepton_mass_slot_from_index(WCoefMapper::lepton_index_from_cq1(coef));
    matching_info[QCDOrder::LO] = {
        {}, // sources seront insérées ensuite par boucle
        [lepton_mass_slot](const ParamSrc& src) { return compute_LO(src, lepton_mass_slot); },
        get_lhaid_from_name(QCDOrder::LO)
    };

    auto& sources = matching_info[QCDOrder::LO].sources;


    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 2});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 3});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 19});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 6});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 7});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 8});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 9});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 11});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_SM", 4});
    // sources.insert({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1});
    sources.insert({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,0}});
    sources.insert({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,1}});
    sources.insert({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,2}});
    sources.insert({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}});
    sources.insert({ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}});
    sources.insert({ParameterType::WILSON, "WPARAM_MATCH_SM", 6});
    // NMSSMScalarMatching::compute compares m_a1 with the matching scale.
    sources.insert({ParameterType::WILSON, "EW_SCALE", 1});
    sources.insert({ParameterType::SM, "SMINPUTS", 2});
    sources.insert({ParameterType::SM, "MASS", 3});
    sources.insert({ParameterType::SM, "MASS", 24});
    sources.insert({ParameterType::SM, "MASS", 25});
    sources.insert({ParameterType::SM, "GAUGE", 2});
    sources.insert({ParameterType::BSM, "ALPHA", LhaID()});
    sources.insert({ParameterType::BSM, "MASS", 37});
    sources.insert({ParameterType::BSM, "MASS", 45});
    sources.insert({ParameterType::BSM, "MASS", 46});
    sources.insert({ParameterType::BSM, "MASS", 2000013});
    sources.insert({ParameterType::BSM, "MASS", 1000006});
    sources.insert({ParameterType::BSM, "MASS", 2000006});
    sources.insert({ParameterType::BSM, "MASS", 25});
    sources.insert({ParameterType::BSM, "MASS", 35});
    sources.insert({ParameterType::BSM, "MASS", 36});
    sources.insert({ParameterType::BSM, "HMIX", 1});
    sources.insert({ParameterType::BSM, "MINPAR", 3});
    sources.insert({ParameterType::WILSON, "EPSILON_SUSY", 5});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {14,0}});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {14,1}});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {14,2}});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {14,3}});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {14,4}});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {14,5}});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {16,0}});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {16,1}});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {16,2}});

    // Boucles dynamiques
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 3; ++j) {
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {1, i, j}});
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {2, i, j}});
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {9, i, j}});
        }
    }

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 6; ++j) {
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {3, i, j, 1}});
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {4, i, j, 2}});
        }
    }

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                sources.insert({ParameterType::WILSON, "MATRIX_BSM", {12, i, 0, j, k}});
                sources.insert({ParameterType::WILSON, "MATRIX_BSM", {12, i, 1, j, k}});
            }
        }
    }

    // UMIX / VMIX with index pairs (2,0), (2,1), ..., (2,9)
    for (int i = 1; i < 3; ++i) {
        sources.insert({ParameterType::BSM, "UMIX", {1, i}});
        sources.insert({ParameterType::BSM, "VMIX", {1, i}});
        sources.insert({ParameterType::BSM, "UMIX", {2, i}});
        sources.insert({ParameterType::BSM, "VMIX", {2, i}});
    }

    for (int i = 0; i < 3; ++i) {
        sources.insert({ParameterType::WILSON, "MATRIX_BSM", {6, 0, i, 1}});
        sources.insert({ParameterType::WILSON, "MATRIX_BSM", {6, 1, i, 1}});
        sources.insert({ParameterType::WILSON, "MATRIX_BSM", {5, 0, i, 1}});
        sources.insert({ParameterType::WILSON, "MATRIX_BSM", {5, 1, i, 1}});
    }

    // HMIX & AMIX matrix elements used in the else branch
    sources.insert({ParameterType::BSM, "NMHMIX", {0+1, 0+1}});
    sources.insert({ParameterType::BSM, "NMHMIX", {0+1, 1+1}});
    sources.insert({ParameterType::BSM, "NMHMIX", {1+1, 0+1}});
    sources.insert({ParameterType::BSM, "NMHMIX", {1+1, 1+1}});
    sources.insert({ParameterType::BSM, "NMHMIX", {2+1, 0+1}});
    sources.insert({ParameterType::BSM, "NMHMIX", {2+1, 1+1}});
            // {ParameterType::BSM, "NMHMIX", {3+1, 0+1}},
    sources.insert({ParameterType::BSM, "NMAMIX" , {0+1, 0+1}});
    sources.insert({ParameterType::BSM, "NMAMIX" , {0+1, 1+1}});
    sources.insert({ParameterType::BSM, "NMAMIX" , {1+1, 0+1}});
    sources.insert({ParameterType::BSM, "NMAMIX" , {1+1, 1+1}});

    matching_info[QCDOrder::NLO] = {
        {}, // les sources sont ajoutées juste après
        [lepton_mass_slot](const ParamSrc& src) { return compute_NLO(src, lepton_mass_slot); },
        get_lhaid_from_name(QCDOrder::NLO)
    };

    auto& sources_NLO = matching_info[QCDOrder::NLO].sources;

    // Insertion directe
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 1});
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 2});
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 3});
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 19});
    // sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 6});
    // sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 7});
    // sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 8});
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 9});
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 11});
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot});
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_SM", 4});
    // sources_NLO.insert({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1});
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}});
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}});
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_MATCH_SM", 6});
    sources_NLO.insert({ParameterType::WILSON, "EW_SCALE", 1});
    // sources_NLO.insert({ParameterType::BSM, "ALPHA", LhaID()});
    sources_NLO.insert({ParameterType::SM, "MASS", 3});
    sources_NLO.insert({ParameterType::SM, "MASS", 24});
    // sources_NLO.insert({ParameterType::SM, "MASS", 25});
    sources_NLO.insert({ParameterType::SM, "SMINPUTS", 4}); //MASS Z in SMINPUTS
    // sources_NLO.insert({ParameterType::SM, "MASS", 23});
    sources_NLO.insert({ParameterType::SM, "GAUGE", 2});
    sources_NLO.insert({ParameterType::SM, "VCKM", {0, 1}});
    sources_NLO.insert({ParameterType::SM, "VCKM", {1, 1}});
    sources_NLO.insert({ParameterType::SM, "VCKM", {2, 1}});
    sources_NLO.insert({ParameterType::SM, "VCKM", {0, 2}});
    sources_NLO.insert({ParameterType::SM, "VCKM", {1, 2}});
    sources_NLO.insert({ParameterType::SM, "VCKM", {2, 2}});
    sources_NLO.insert({ParameterType::BSM, "MASS", 37});
    // sources_NLO.insert({ParameterType::BSM, "MASS", 35});
    sources_NLO.insert({ParameterType::BSM, "MASS", 45});
    sources_NLO.insert({ParameterType::BSM, "MASS", 46});
    sources_NLO.insert({ParameterType::BSM, "HMIX", 1});
    sources_NLO.insert({ParameterType::BSM, "HMIX", 2});
    sources_NLO.insert({ParameterType::WILSON, "EPSILON_SUSY", 5});

    // Boucles : WPARAM_SI_BSM et mixings
    for (int i = 0; i < 6; ++i) {
        
        sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {14, i}});
        

    }

    for (int i = 1; i < 3; ++i) {
        

        sources_NLO.insert({ParameterType::BSM, "UMIX", {1, i}});
        sources_NLO.insert({ParameterType::BSM, "VMIX", {1, i}});
        sources_NLO.insert({ParameterType::BSM, "UMIX", {2, i}});
        sources_NLO.insert({ParameterType::BSM, "VMIX", {2, i}});
        

    }

    // MATRIX_BSM blocs 1, 2, 9
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 3; ++j) {
            sources_NLO.insert({ParameterType::WILSON, "MATRIX_BSM", {1, i, j}});
            sources_NLO.insert({ParameterType::WILSON, "MATRIX_BSM", {2, i, j}});
            
        }
        for (int j = 0; j < 6; ++j) {
            sources_NLO.insert({ParameterType::WILSON, "MATRIX_BSM", {9, i, j}});
        }
    }

    // MATRIX_BSM blocs 3, 4
    for (int i = 0; i < 2; ++i) {
        sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {13, i}});
        for (int j = 0; j < 6; ++j) {
            sources_NLO.insert({ParameterType::WILSON, "MATRIX_BSM", {3, i, j, 1}});
            sources_NLO.insert({ParameterType::WILSON, "MATRIX_BSM", {4, i, j, 2}});
            sources_NLO.insert({ParameterType::WILSON, "MATRIX_BSM", {5, i, j, 1}});
            sources_NLO.insert({ParameterType::WILSON, "MATRIX_BSM", {5, i, j, 2}});
            sources_NLO.insert({ParameterType::WILSON, "MATRIX_BSM", {6, i, j, 1}});
            sources_NLO.insert({ParameterType::WILSON, "MATRIX_BSM", {6, i, j, 2}});
        }
    }

    // MATRIX_BSM bloc 12
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                sources_NLO.insert({ParameterType::WILSON, "MATRIX_BSM", {12, i, 0, j, k}});
                sources_NLO.insert({ParameterType::WILSON, "MATRIX_BSM", {12, i, 1, j, k}});
                // sources_NLO.insert({ParameterType::WILSON, "MATRIX_BSM", {12, i, 2, j, k}});
            }
        }
    }

    // // RECKM i=1→3, j=1→3 → 11,12,13,21,22,23,31,32,33
    // for (int i = 1; i <= 3; ++i) {
    //     for (int j = 1; j <= 3; ++j) {
    //         sources_NLO.insert({ParameterType::SM, "RECKM", i * 10 + j});
    //     }
    // }

    // WPARAM_MATCH_BSM indexé
    for (int fe = 0; fe < 3; ++fe) {
        sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {16, fe}});
        sources_NLO.insert({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe}});
    }

    // HMIX & AMIX matrix elements used in the else branch
    sources_NLO.insert({ParameterType::BSM, "NMHMIX", {0+1, 0+1}});
    sources_NLO.insert({ParameterType::BSM, "NMHMIX", {0+1, 1+1}});
    sources_NLO.insert({ParameterType::BSM, "NMHMIX", {1+1, 0+1}});
    sources_NLO.insert({ParameterType::BSM, "NMHMIX", {1+1, 1+1}});
    sources_NLO.insert({ParameterType::BSM, "NMHMIX", {2+1, 0+1}});
    sources_NLO.insert({ParameterType::BSM, "NMHMIX", {2+1, 1+1}});
            // {ParameterType::BSM, "NMHMIX", {3+1, 0+1}},
    sources_NLO.insert({ParameterType::BSM, "NMAMIX" , {0+1, 0+1}});
    sources_NLO.insert({ParameterType::BSM, "NMAMIX" , {0+1, 1+1}});
    sources_NLO.insert({ParameterType::BSM, "NMAMIX" , {1+1, 0+1}});
    sources_NLO.insert({ParameterType::BSM, "NMAMIX" , {1+1, 1+1}});

}

scalar_t CQ1_susy::compute_LO(const ParamSrc& src) {
    return compute_LO(src, WCoefMapper::lepton_mass_slot_from_index(1));
}

scalar_t CQ1_susy::compute_LO(const ParamSrc& src, int lepton_mass_slot) {
    const auto nmssm = nmssm_scalar_matching::compute(src, lepton_mass_slot);
    if (nmssm.active) return nmssm.cq1;

    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1});

    double lu = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 7);
    double ld = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 8);
    double aY = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 11);
    // double yt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    // double Q_match = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);;
    double mW = src.get_val(ParameterType::SM, "MASS", 24);
    double mh0 = src.get_val(ParameterType::SM, "MASS", 25);
    double mH = src.get_val(ParameterType::BSM, "MASS", 37);
    double mH0 = src.get_val(ParameterType::BSM, "MASS", 35);
    double sw2 = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);
    double mass_b_muW = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1});
    double g2 = src.get_val(ParameterType::SM, "GAUGE", 2);
    double tanb = src.get_val(ParameterType::BSM, "MINPAR", 3);

    double muQ = src.get_val(ParameterType::BSM, "HMIX", 1);
    double xh = pow(mh0/mW, 2);
    double xH = pow(mH/mW, 2);
    double xH0 = pow(mH0/mW, 2);
    double alpha = src.get_val(ParameterType::BSM, "ALPHA", LhaID());
    double beta = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 6);
    double kappa = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 19);
    double epsfac = src.get_val(ParameterType::WILSON, "EPSILON_SUSY", 5);
    complex_t BQ10c1=0.;
    complex_t BQ10c2=0.;

    double Dp, Dm;
    // double a0a{0}, a0b{0}, a0c, a0Q1{0}, a0Q2{0};
    double a0a{0}, a0b{0}, a0c, a0Q1{0};
    double a1{0};

    complex_t NQ10c=0.;
    for(int ie=0;ie<2;ie++) {
        for(int je=0;je<2;je++) {
            for(int ae=0;ae<6;ae++) {
                for(int be=0;be<3;be++) { 
                    BQ10c1+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, je, ae, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6, ie, be, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5, je, be, 1})*f50(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)));
                    BQ10c2+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, je, ae, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5, ie, be, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6, je, be, 1})*fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}))*f60(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)));
                    for(int me=0;me<6;me++) {
                        for(int ne=0;ne<3;ne++) {
                            Dp=0.;
                            Dm=0.;
                            for(int fe=0;fe<3;fe++) 
                            {
                                Dp+=src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe})/sqrt(2.)/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*muQ*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, me, fe})+src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, me, fe}));
                                Dm+=src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe})/sqrt(2.)/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*muQ*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, me, fe})-src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, me, fe}));
                            }
                            a0a=-(fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}))*f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.))*src.get_val(ParameterType::BSM, "UMIX", {ie+1,1+1})*src.get_val(ParameterType::BSM, "VMIX", {je+1, 0+1}))*kron(ae,me);
                            a0b=-(f40(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.))*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*src.get_val(ParameterType::BSM, "VMIX", {ie+1, 0+1}))*kron(ae,me);
                            a0c=1./mW*f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, me})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))*kron(ie,je);
                            a0Q1=a0a+a0b+Dp*a0c;
                            
                            a1=src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/sqrt(2.)/mW*f80(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))*kron(ie,je)*kron(ae,me);
                            NQ10c+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, me, be})*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*(a0Q1+a1*tanb);
                        }
                    
                    }
                }
            }
        }
    }
    complex_t BQ10c=(BQ10c1+BQ10c2)*kappa*mW*mW/2./g2/g2/sw2;
    // printf("BQ10c1 : %.14lf\n", BQ10c1);
    // printf("BQ10c2 : %.14lf\n", BQ10c2);
    NQ10c*=src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot)*(tanb)*tanb/mW/(mH*mH-mW*mW)*aY*mass_b_muW/sw2;
    double le = -tanb;
    double G1=-3./4.+ld*lu*F4SP(xt,xH)+lu*lu*F5SP(xt,xH);
    double G2=ld*(ld*lu+1.)*F6SP(xt,xH)-ld*lu*lu*F7SP(xt,xH)
    +lu*lu*(ld*F8SP(xt,xH)+lu*F9SP(xt,xH)-lu*F10SP(xt,xH))+lu*F11SP(xt,xH)-lu*F12SP(xt,xH);

    double CSn_2HDM=xt*(F0SP(xt)+le*(ld*F1SP(xt,xH)+lu*F2SP(xt,xH))+le*lu*F3SP(xt,xH))
    +xt/2./xh*(sin(alpha-beta)+cos(alpha-beta)*le)*(sin(alpha-beta)*G1+cos(alpha-beta)*G2)
    +xt/2./xH0*(cos(alpha-beta)-sin(alpha-beta)*le)*(cos(alpha-beta)*G1-sin(alpha-beta)*G2);
    complex_t CQ1H_0=CSc_2HDM(xH,xt,lu,ld,le)+CSn_2HDM;
    CQ1H_0*=(src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot)*mass_b_muW/mW/mW)/sw2;

    // printf("sw2 = %.9lf\n", sw2);
	// printf("mass_b_muW = %.9lf\n", mass_b_muW);
	// printf("param->mass_W = %.9lf\n", mW);
    // printf("CSn_2HDM : %.9lf\n", CSn_2HDM);
    // printf("xH : %.9lf\n", xH);
    // printf("xt : %.9lf\n", xt);
    // printf("lu : %.9lf\n", lu);
    // printf("ld : %.9lf\n", ld);
    // printf("le : %.9lf\n", le);
    // printf("alpha : %.9lf\n", alpha);
    // printf("beta : %.9lf\n", beta);
    // printf("G1 : %.9lf\n", G1);
    // printf("G2 : %.9lf\n", G2);
    // printf("CSc_2HDM(xH,xt,lu,ld,le) : %.9lf\n", CSc_2HDM(xH,xt,lu,ld,le));
    
    // printf("NQ10c : %.9lf\n", NQ10c.real());
	// printf("BQ10c : %.9lf\n", BQ10c.real());
    complex_t CQ1charg_0=NQ10c+BQ10c;
    complex_t coeff_temp = (CQ1charg_0+CQ1H_0)/epsfac;
    // printf("kappa : %.9lf\n", kappa);
    // std::cout << "CQ1charg_0 " << CQ1charg_0 << std::endl;
    // std::cout << "CQ1H_0 " << CQ1H_0.real() << std::endl;
    return coeff_temp;
}

scalar_t CQ1_susy::compute_NLO(const ParamSrc& src) {
    return compute_NLO(src, WCoefMapper::lepton_mass_slot_from_index(1));
}

scalar_t CQ1_susy::compute_NLO(const ParamSrc& src, int lepton_mass_slot) {
    if (nmssm_scalar_matching::is_active()) return scalar_t(0.0);

    // LOG_INFO("CQ1_susy::NLO 1");
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1});

    // double lu = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 7);
    // double ld = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 8);
    double z = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 1);
    double aY = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 11);
    // double yt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    double Q_match = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);;
    double mW = src.get_val(ParameterType::SM, "MASS", 24);
    // double mh0 = src.get_val(ParameterType::SM, "MASS", 25);
    double mH = src.get_val(ParameterType::BSM, "MASS", 37);
    // double mH0 = src.get_val(ParameterType::BSM, "MASS", 35);
    double sw2 = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);
    double mass_b_muW = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1});
    double mass_top_muW = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 6);
    double g2 = src.get_val(ParameterType::SM, "GAUGE", 2);
    double tanb = src.get_val(ParameterType::BSM, "HMIX", 2);

    double muQ = src.get_val(ParameterType::BSM, "HMIX", 1);
    // double xh = pow(mh0/mW, 2);
    // double xH = pow(mH/mW, 2);
    // double xH0 = pow(mH0/mW, 2);
    // double alpha = src.get_val(ParameterType::BSM, "ALPHA", LhaID());
    // double beta = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 6);

    double kappa = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 19);
    double epsfac = src.get_val(ParameterType::WILSON, "EPSILON_SUSY", 5);

    complex_t NQ11H=-src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot)*(tanb)*tanb/4./mW/mW*(f141(xt,z)+8.*xt*(f30(xt,z)+xt*(f30(xt*1.0001,z)-f30(xt*0.9999,z))/0.0002)*log(Q_match*Q_match/mass_top_muW/mass_top_muW));
    complex_t BQ11H=src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot)*(tanb)*tanb/4./mW/mW*(f111(xt,z)+8.*(f70(xt*1.0001,z)-f70(xt*0.9999,z))/0.0002*log(Q_match*Q_match/mass_top_muW/mass_top_muW));
    complex_t CQ1H_1=(NQ11H+BQ11H)*mass_b_muW/sw2;
    // complex_t CQ2H_1=-CQ1H_1;

    // LOG_INFO("CQ1_susy::NLO 2");
    complex_t BQ11c1=0.;
    complex_t BQ11c2=0.;
    complex_t NQ11c=0.;
    complex_t NQ21c=0.;
    complex_t BQ11f1=0.;
    complex_t BQ11f2=0.;
    complex_t NQ11f=0.;
    complex_t NQ21f=0.;

    complex_t Dp{0}, Dm{0}, temp{0}, temp2{0};
    complex_t a0a{0}, a0b{0}, a0c{0}, a0Q1{0}, a0Q2{0}, a0p{0}, a1{0},a2p{0};

    for(int ie=0;ie<2;ie++) {
        for(int ae=0;ae<6;ae++){
            for(int je=0;je<2;je++)  {
                for(int be=0;be<3;be++) {
                    BQ11c1+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, je, ae, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6, ie, be, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5, je, be, 1})*(f121(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))+4.*(f50(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*1.0001,pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))-f50(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*0.9999,pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)))/0.0002*log(pow(mass_top_muW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}),2.))));
                    BQ11c2+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, je, ae, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5, ie, be, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6, je, be, 1})*fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}))*(f131(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))+4.*(f60(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*1.0001,pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))-f60(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*0.9999,pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)))/0.0002*log(pow(mass_top_muW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}),2.))));
                    // LOG_INFO("CQ1_susy::NLO 2.-1");
                    for(int me=0;me<6;me++){ 
                        for(int ne=0;ne<3;ne++) {
                            Dp=0.;
                            Dm=0.;
                            for(int fe=0;fe<3;fe++) { 	
                                Dp+=src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe})/sqrt(2.)/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*muQ*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, me, fe})+src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, me, fe}));
                                Dm+=src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe})/sqrt(2.)/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*muQ*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, me, fe})-src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, me, fe}));
                            }
                            a0a=-(fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}))*(f181(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.))+4.*(f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.)*1.0001)-f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.)*0.9999))/0.0002*log(pow(mass_top_muW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}),2.)))*src.get_val(ParameterType::BSM, "UMIX", {ie+1, 1+1})*src.get_val(ParameterType::BSM, "VMIX", {je+1, 0+1}))*kron(ae,me);
                            a0b=-((f191(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.))+4.*(f40(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.)*1.0001)-f40(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.)*0.9999))/0.0002*log(pow(mass_top_muW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}),2.)))*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*src.get_val(ParameterType::BSM, "UMIX", {ie+1, 0+1}))*kron(ae,me);
                            a0c=1./mW*(f171(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, me})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))+4.*(f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, me})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))+(f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*1.0001,pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, me})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))-f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*0.9999,pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, me})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)))/0.0002+(f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, me})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*1.0001)-f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, me})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*0.9999))/0.0002)*log(pow(mass_top_muW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}),2.)))*kron(ie,je);
                        
                            a0Q1=a0a+a0b+Dp*a0c;
                            a0Q2=-a0a+a0b+Dm*a0c;
                            a0p=4.*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne})/mW/(src.get_val(ParameterType::SM, "VCKM", {be, 2})*src.get_val(ParameterType::SM, "VCKM", {ne, 1})/src.get_val(ParameterType::SM, "VCKM", {2, 2})/src.get_val(ParameterType::SM, "VCKM", {2, 1}))/src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*f151(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))*kron(ie,je)*kron(ae,me)*kron(be,ne);
                            a1=src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/sqrt(2.)/mW*(f161(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))+4.*(f80(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*1.0001)-f80(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*0.9999))/0.0002*log(pow(mass_top_muW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}),2.)))*kron(ie,je)*kron(ae,me);
                            a2p=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, me, be})*(src.get_val(ParameterType::SM, "VCKM", {be, 2})*src.get_val(ParameterType::SM, "VCKM", {ne, 1})/src.get_val(ParameterType::SM, "VCKM", {2, 2})/src.get_val(ParameterType::SM, "VCKM", {2, 1}))*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})/2./mW*f151(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))*kron(ie,je)*kron(ae,me)*kron(be,ne);
                            
                            NQ11c+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, me, be})*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*(a0Q1+a1*tanb)
                            +src.get_val(ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne})*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*a0p
                            +src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, me, be})*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*a2p*pow(src.get_val(ParameterType::SM, "MASS", 3)*tanb,2.);	
                            NQ21c+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, me, be})*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*(a0Q2+a1*tanb)
                            +src.get_val(ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne})*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*a0p
                            +src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, me, be})*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*a2p*pow(src.get_val(ParameterType::SM, "MASS", 3)*tanb,2.);
                        }
                    }
                
                }
                // LOG_INFO("CQ1_susy::NLO 2.0");
                for(int be=0;be<6;be++) {
                    for(int ce=0;ce<6;ce++) {
                        for(int fe=0;fe<3;fe++) {
                            BQ11f1+=-src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, je, be, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2})*pow(mW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9, ae, ce})*src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9, ce, be})*(1.+log(pow(mass_top_muW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}),2.)))	*(f90(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, fe})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6, ie, fe, 2})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5, je, fe, 2}));
                            BQ11f2+=-src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, je, be, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2})*pow(mW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9, ae, ce})*src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9, ce, be})*(1.+log(pow(mass_top_muW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}),2.)))	*(fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}))*f100(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, fe})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5, ie, fe, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6, je, fe, 1}));

                        }
                    }
                }
            }

            // LOG_INFO("CQ1_susy::NLO 2.1");

            for(int me=0;me<3;me++) {
                for(int ne=0;ne<3;ne++) {
                    for(int de=0;de<6;de++) {
                        for(int ke=0;ke<6;ke++) {
                            
                            temp2 = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, me, ne})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, de, me})*src.get_val(ParameterType::BSM, "UMIX", {ie+1, 1+1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9, ae, ke})*src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9, ke, de})*(1.+log(pow(mass_top_muW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}),2.)))*tanb*src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/sqrt(2.)*f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.));
                            NQ11f+=temp2;
                            NQ21f+=temp2;
                            for(int ce=0;ce<6;ce++) {
                                Dp=0.;
                                Dm=0.;
                                for(int fe=0;fe<3;fe++) 
                                {		
                                    Dp+=src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe})/sqrt(2.)/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*muQ*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, de, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ce, fe})+src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, de, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, ce, fe})); 
                                    Dm+=src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe})/sqrt(2.)/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*muQ*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, de, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ce, fe})-src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, de, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, ce, fe})); 
                                }
                                temp=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, me, ne})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ce, me})*src.get_val(ParameterType::BSM, "UMIX", {ie+1, 1+1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9, ae, ke})*src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9, ke, de})*
                                (1.+log(pow(mass_top_muW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}),2.)))*f60(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.));	
                                NQ11f+=Dp*temp;
                                NQ21f+=Dm*temp;
                            }
                            for(int je=0;je<2;je++) {
                                temp=-src.get_val(ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, me, ne})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, de, me})*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9, ae, ke})*src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9, ke, de})*
                                (1.+log(pow(mass_top_muW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}),2.)))*mW*(fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}))*f60(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.))*src.get_val(ParameterType::BSM, "UMIX", {ie+1, 1+1})*src.get_val(ParameterType::BSM, "VMIX", {je+1, 0+1})+
                                f50(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.))*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*src.get_val(ParameterType::BSM, "VMIX", {ie+1, 0+1})); 
                                NQ11f+=temp;
                                NQ21f+=-temp;
                            }

                        }
                    }
                }
                    
            }
        }
    }

    // LOG_INFO("CQ1_susy::NLO 3");
    
    complex_t BQ11c=(BQ11c1+BQ11c2)*kappa*mW*mW/2./g2/g2/sw2;
    // complex_t BQ21c=-(BQ11c1-BQ11c2)*kappa*mW*mW/2./g2/g2/sw2;


    NQ11c*=src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot)*(tanb)*tanb/mW/(mH*mH-mW*mW)*aY*mass_b_muW/sw2;
    NQ21c*=-src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot)*(tanb)*tanb/mW/(mH*mH-mW*mW)*aY*mass_b_muW/sw2;
    
    complex_t CQ1charg_1=NQ11c+BQ11c;
    


    
    // complex_t CQ2charg_1=NQ21c+BQ21c;

    complex_t BQ11f=(BQ11f1+BQ11f2)*2./3.*kappa/g2/g2/sw2;
    // complex_t BQ21f=-(BQ11f1-BQ11f2)*2./3.*kappa/g2/g2/sw2;
    
    
    NQ11f*=-4./3.*src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot)*(tanb)*tanb/mW/mW/(mH*mH-mW*mW)*aY*mass_b_muW/sw2;

    NQ21f*=4./3.*src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot)*(tanb)*tanb/mW/mW/(mH*mH-mW*mW)*aY*mass_b_muW/sw2;


    complex_t CQ1four_1=NQ11f+BQ11f;
    complex_t coeff_temp = (CQ1H_1+CQ1charg_1)/epsfac+CQ1four_1;

    // LOG_INFO("CQ1_susy::NLO 4");


    return coeff_temp;
}

CQ2_susy::CQ2_susy(WCoef coef) : WilsonCoefficient(WCoefMapper::str(coef) + "_SUSY", GroupMapper::str(WGroup::BScalar, ScaleType::MATCHING)) {
    const int lepton_mass_slot = WCoefMapper::lepton_mass_slot_from_index(WCoefMapper::lepton_index_from_cq2(coef));
    matching_info[QCDOrder::LO] = {
        {}, // sources
        [lepton_mass_slot](const ParamSrc& src) { return compute_LO(src, lepton_mass_slot); },
        get_lhaid_from_name(QCDOrder::LO)
    };

    auto& sources = matching_info[QCDOrder::LO].sources;

    // Insertion directe
    // sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 1});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 2});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 3});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 4});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 19});
    // sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 6});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 7});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 8});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 9});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 11});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_SM", 4});
    sources.insert({ParameterType::WILSON, "WPARAM_SI_SM", 5});
    sources.insert({ParameterType::WILSON, "WPARAM_RUN_SM", 1});
    // sources.insert({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1});
    sources.insert({ParameterType::WILSON, "EW_SCALE", 1});
    // sources.insert({ParameterType::BSM, "ALPHA", LhaID()});
    sources.insert({ParameterType::BSM, "MASS", 35});
    sources.insert({ParameterType::SM, "MASS", 25});
    sources.insert({ParameterType::BSM, "MASS", 37});
    sources.insert({ParameterType::BSM, "MASS", 45});
    sources.insert({ParameterType::BSM, "MASS", 46});
    sources.insert({ParameterType::BSM, "MASS", 2000013});
    sources.insert({ParameterType::BSM, "MASS", 1000006});
    sources.insert({ParameterType::BSM, "MASS", 2000006});
    // sources.insert({ParameterType::BSM, "MASS", 25});
    // sources.insert({ParameterType::BSM, "MASS", 35});
    sources.insert({ParameterType::BSM, "MASS", 36});
    sources.insert({ParameterType::BSM, "HMIX", 1});
    sources.insert({ParameterType::BSM, "MINPAR", 3});
    sources.insert({ParameterType::WILSON, "EPSILON_SUSY", 5});
    // sources.insert({ParameterType::SM, "MASS", 23});
    sources.insert({ParameterType::SM, "SMINPUTS", 4}); //MASS Z in SMINPUTS
    sources.insert({ParameterType::SM, "MASS", 24});
    sources.insert({ParameterType::SM, "GAUGE", 2});
    sources.insert({ParameterType::SM, "QCD", LhaID(5,1)});
    sources.insert({ParameterType::SM, "SMINPUTS", 2});

    // MATRIX_BSM 3 & 4
    for (int je = 0; je < 2; ++je) {
        for (int ae = 0; ae < 6; ++ae) {
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {3, je, ae, 1}});
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {4, je, ae, 2}});
        }
    }

    // MATRIX_BSM 5 & 6
    for (int ie = 0; ie < 2; ++ie) {
        for (int be = 0; be < 3; ++be) {
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {5, ie, be, 1}});
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {6, ie, be, 1}});
        }
    }

    // WPARAM_SI_BSM mixings
    for (int ie = 0; ie < 2; ++ie) {
        sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}});
    }
    for (int ae = 0; ae < 6; ++ae) {
        sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}});
    }
    for (int be = 0; be < 3; ++be) {
        sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}});
    }

    for (int i = 1; i < 3; ++i) {
        sources.insert({ParameterType::BSM, "UMIX", {1, i}});
        sources.insert({ParameterType::BSM, "VMIX", {1, i}});
        sources.insert({ParameterType::BSM, "UMIX", {2, i}});
        sources.insert({ParameterType::BSM, "VMIX", {2, i}});

    }

    // WPARAM_MATCH_SM
    sources.insert({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}});
    sources.insert({ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}});
    sources.insert({ParameterType::WILSON, "WPARAM_MATCH_SM", 6});

    // WPARAM_MATCH_BSM indexés
    for (int fe = 0; fe < 3; ++fe) {
        sources.insert({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe}});
    }

    // MATRIX_BSM 1 & 2
    for (int ae = 0; ae < 6; ++ae) {
        for (int fe = 0; fe < 3; ++fe) {
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe}});
        }
    }
    for (int me = 0; me < 6; ++me) {
        for (int be = 0; be < 3; ++be) {
            sources.insert({ParameterType::WILSON, "MATRIX_BSM", {1, me, be}});
        }
    }

    // MATRIX_BSM 12
    for (int ae = 0; ae < 6; ++ae) {
        for (int ie = 0; ie < 2; ++ie) {
            for (int be = 0; be < 3; ++be) {
                for (int ne = 0; ne < 3; ++ne) {
                    sources.insert({ParameterType::WILSON, "MATRIX_BSM", {12, ae, ie, be, ne}});
                }
            }
        }
    }

    // HMIX & AMIX matrix elements used in the else branch
    sources.insert({ParameterType::BSM, "NMHMIX", {0+1, 0+1}});
    sources.insert({ParameterType::BSM, "NMHMIX", {0+1, 1+1}});
    sources.insert({ParameterType::BSM, "NMHMIX", {1+1, 0+1}});
    sources.insert({ParameterType::BSM, "NMHMIX", {1+1, 1+1}});
    sources.insert({ParameterType::BSM, "NMHMIX", {2+1, 0+1}});
    sources.insert({ParameterType::BSM, "NMHMIX", {2+1, 1+1}});
            // {ParameterType::BSM, "NMHMIX", {3+1, 0+1}},
    sources.insert({ParameterType::BSM, "NMAMIX" , {0+1, 0+1}});
    sources.insert({ParameterType::BSM, "NMAMIX" , {0+1, 1+1}});
    sources.insert({ParameterType::BSM, "NMAMIX" , {1+1, 0+1}});
    sources.insert({ParameterType::BSM, "NMAMIX" , {1+1, 1+1}});

    matching_info[QCDOrder::NLO] = {
        {}, // sources
        [lepton_mass_slot](const ParamSrc& src) { return compute_NLO(src, lepton_mass_slot); },
        get_lhaid_from_name(QCDOrder::NLO)
    };

    auto& sources_NLO = matching_info[QCDOrder::NLO].sources;

    // Insertion directe
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 1});
    // sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 2});
    // sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 3});
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 19});
    // sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 6});
    // sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 7});
    // sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 8});
    // sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 9});
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", 11});
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot});
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_SM", 4});
    // sources_NLO.insert({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1});
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}});
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}});
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_MATCH_SM", 6});
    sources_NLO.insert({ParameterType::WILSON, "EW_SCALE", 1});
    sources_NLO.insert({ParameterType::WILSON, "EPSILON_SUSY", 5});
    sources_NLO.insert({ParameterType::SM, "MASS", 3});
    sources_NLO.insert({ParameterType::SM, "MASS", 24});
    sources_NLO.insert({ParameterType::SM, "GAUGE", 2});
    sources_NLO.insert({ParameterType::SM, "VCKM", {0, 1}});
    sources_NLO.insert({ParameterType::SM, "VCKM", {1, 1}});
    sources_NLO.insert({ParameterType::SM, "VCKM", {2, 1}});
    sources_NLO.insert({ParameterType::SM, "VCKM", {0, 2}});
    sources_NLO.insert({ParameterType::SM, "VCKM", {1, 2}});
    sources_NLO.insert({ParameterType::SM, "VCKM", {2, 2}});
    sources_NLO.insert({ParameterType::BSM, "MASS", 37});
    sources_NLO.insert({ParameterType::BSM, "MASS", 45});
    sources_NLO.insert({ParameterType::BSM, "MASS", 46});
    sources_NLO.insert({ParameterType::BSM, "HMIX", 1});
    sources_NLO.insert({ParameterType::BSM, "MINPAR", 3});

    // WPARAM_SI_BSM
    for (int ie = 0; ie < 2; ++ie) {
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}});
    }
    for (int ae = 0; ae < 6; ++ae) {
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}});
    }
    for (int be = 0; be < 3; ++be) {
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}});
    }

    for (int i = 1; i < 3; ++i) {
        sources_NLO.insert({ParameterType::BSM, "UMIX", {1, i}});
        sources_NLO.insert({ParameterType::BSM, "VMIX", {1, i}});
        sources_NLO.insert({ParameterType::BSM, "UMIX", {2, i}});
        sources_NLO.insert({ParameterType::BSM, "VMIX", {2, i}});

    }

    // MATRIX_BSM: blocs 3 & 4
    for (int je = 0; je < 2; ++je) {
    for (int ae = 0; ae < 6; ++ae) {
        sources_NLO.insert({ParameterType::WILSON, "MATRIX_BSM", {3, je, ae, 1}});
        sources_NLO.insert({ParameterType::WILSON, "MATRIX_BSM", {4, je, ae, 2}});
    }
    }

    // MATRIX_BSM: blocs 5 & 6
    for (int ie = 0; ie < 2; ++ie) {
    for (int be = 0; be < 3; ++be) {
        sources_NLO.insert({ParameterType::WILSON, "MATRIX_BSM", {5, ie, be, 1}});
        sources_NLO.insert({ParameterType::WILSON, "MATRIX_BSM", {5, ie, be, 2}});
        sources_NLO.insert({ParameterType::WILSON, "MATRIX_BSM", {6, ie, be, 1}});
        sources_NLO.insert({ParameterType::WILSON, "MATRIX_BSM", {6, ie, be, 2}});
    }
    }

    // MATRIX_BSM: bloc 2 et 1
    for (int ae = 0; ae < 6; ++ae) {
    for (int fe = 0; fe < 3; ++fe) {
        sources_NLO.insert({ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe}});
        sources_NLO.insert({ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe}});
    }
    }

    // MATRIX_BSM: bloc 12
    for (int ae = 0; ae < 6; ++ae) {
    for (int ie = 0; ie < 2; ++ie) {
        for (int me = 0; me < 3; ++me) {
            for (int ne = 0; ne < 3; ++ne) {
                sources_NLO.insert({ParameterType::WILSON, "MATRIX_BSM", {12, ae, ie, me, ne}});
            }
        }
    }
    }

    // MATRIX_BSM: bloc 9
    for (int ae = 0; ae < 6; ++ae) {
    for (int ce = 0; ce < 6; ++ce) {
        sources_NLO.insert({ParameterType::WILSON, "MATRIX_BSM", {9, ae, ce}});
    }
    }

    // WPARAM_MATCH_BSM indices {2, fe}
    for (int fe = 0; fe < 3; ++fe) {
    sources_NLO.insert({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe}});
    }

    // HMIX & AMIX matrix elements used in the else branch
    sources_NLO.insert({ParameterType::BSM, "NMHMIX", {0+1, 0+1}});
    sources_NLO.insert({ParameterType::BSM, "NMHMIX", {0+1, 1+1}});
    sources_NLO.insert({ParameterType::BSM, "NMHMIX", {1+1, 0+1}});
    sources_NLO.insert({ParameterType::BSM, "NMHMIX", {1+1, 1+1}});
    sources_NLO.insert({ParameterType::BSM, "NMHMIX", {2+1, 0+1}});
    sources_NLO.insert({ParameterType::BSM, "NMHMIX", {2+1, 1+1}});
            // {ParameterType::BSM, "NMHMIX", {3+1, 0+1}},
    sources_NLO.insert({ParameterType::BSM, "NMAMIX" , {0+1, 0+1}});
    sources_NLO.insert({ParameterType::BSM, "NMAMIX" , {0+1, 1+1}});
    sources_NLO.insert({ParameterType::BSM, "NMAMIX" , {1+1, 0+1}});
    sources_NLO.insert({ParameterType::BSM, "NMAMIX" , {1+1, 1+1}});

    // // RECKM approximatif (0 à 29)
    // for (int i = 0; i < 30; ++i) {
    // sources_NLO.insert({ParameterType::SM, "RECKM", i});
    // }

}

scalar_t CQ2_susy::compute_LO(const ParamSrc& src) {
    return compute_LO(src, WCoefMapper::lepton_mass_slot_from_index(1));
}

scalar_t CQ2_susy::compute_LO(const ParamSrc& src, int lepton_mass_slot) {
    const auto nmssm = nmssm_scalar_matching::compute(src, lepton_mass_slot);
    if (nmssm.active) return nmssm.cq2;

    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1});

    double lu = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 7);
    double ld = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 8);
    // double z = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 1);
    double aY = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 11);
    // double yt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    double mW = src.get_val(ParameterType::SM, "MASS", 24);
    // double mh0 = src.get_val(ParameterType::SM, "MASS", 25);
    double mH = src.get_val(ParameterType::BSM, "MASS", 37);
    // double mH0 = src.get_val(ParameterType::BSM, "MASS", 35);
    double mA0 = src.get_val(ParameterType::BSM, "MASS", 36);
    double sw2 = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);
    double mass_b_muW = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1});
    double g2 = src.get_val(ParameterType::SM, "GAUGE", 2);
    double tanb = src.get_val(ParameterType::BSM, "MINPAR", 3);

    double muQ = src.get_val(ParameterType::BSM, "HMIX", 1);
    // double xh = pow(mh0/mW, 2);
    double xH = pow(mH/mW, 2);
    // double xH0 = pow(mH0/mW, 2);
    double xA = pow(mA0/mW, 2);
    // double alpha = src.get_val(ParameterType::BSM, "ALPHA", LhaID());
    // double beta = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 6);

    double kappa = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 19);
    double epsfac = src.get_val(ParameterType::WILSON, "EPSILON_SUSY", 5);

    
    complex_t BQ10c1=0.;
    complex_t BQ10c2=0.;

    complex_t Dp, Dm;
    complex_t a0a{0}, a0b{0}, a0c, a0Q1{0}, a0Q2{0};
    complex_t a1{0};
    complex_t NQ20c=0.;

    for(int ie=0;ie<2;ie++) {
        for(int je=0;je<2;je++) {
            for(int ae=0;ae<6;ae++) {
                for(int be=0;be<3;be++) { 
                    BQ10c1+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, je, ae, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6,ie, be, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5, je, be, 1})*f50(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)));
                    BQ10c2+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, je, ae, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5, ie, be, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6, je, be, 1})*fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}))*f60(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)));

                    for(int me=0;me<6;me++) {
                        for(int ne=0;ne<3;ne++) {
                            Dp=0.;
                            Dm=0.;
                            for(int fe=0;fe<3;fe++) 
                            {
                                Dp+=src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe})/sqrt(2.)/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*muQ*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, me, fe})+src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, me, fe}));
                                Dm+=src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe})/sqrt(2.)/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*muQ*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, me, fe})-src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, me, fe}));
                            }
                            a0a=-(fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}))*f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.))*src.get_val(ParameterType::BSM, "UMIX", {ie+1, 1+1})*src.get_val(ParameterType::BSM, "VMIX", {je+1, 0+1}))*kron(ae,me);
                            a0b=-(f40(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.))*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*src.get_val(ParameterType::BSM, "VMIX", {ie+1, 0+1}))*kron(ae,me);
                            a0c=1./mW*f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, me})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))*kron(ie,je);
                            a0Q2=-a0a+a0b+Dm*a0c;
                            
                            a1=src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/sqrt(2.)/mW*f80(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))*kron(ie,je)*kron(ae,me);
                            NQ20c+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, me, be})*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*(a0Q2+a1*tanb);
                        }
                    
                    }
                }
            }
        }
    }
    complex_t BQ20c=-(BQ10c1-BQ10c2)*kappa*mW*mW/2./g2/g2/sw2;
    
    double le = -tanb;
    double G3=ld*(ld*lu+1.)*F6SP(xt,xH)+ld*lu*lu*F7SP(xt,xH)
    +lu*lu*(ld*F8SP(xt,xH)+lu*F9SP(xt,xH)+lu*F10SP(xt,xH))+lu*F11SP(xt,xH)+lu*F12SP(xt,xH);
    double CPn_2HDM=xt*(-le*(ld*F1SP(xt,xH)+lu*F2SP(xt,xH))+le*lu*F3SP(xt,xH))+xt/2./xA*(le)*G3;

    double CQ2H_0=CPc_2HDM(xH,xt,lu,ld,le,sw2)+CPn_2HDM;
    CQ2H_0*=(src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot)*mass_b_muW/mW/mW)/sw2;

    NQ20c*=-src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot)*(tanb)*tanb/mW/(mH*mH-mW*mW)*aY*mass_b_muW/sw2;

    complex_t CQ2charg_0=NQ20c+BQ20c;


    complex_t coeff_temp = (CQ2charg_0+CQ2H_0)/epsfac;

    return coeff_temp;
}


scalar_t CQ2_susy::compute_NLO(const ParamSrc& src) {
    return compute_NLO(src, WCoefMapper::lepton_mass_slot_from_index(1));
}

scalar_t CQ2_susy::compute_NLO(const ParamSrc& src, int lepton_mass_slot) {
    if (nmssm_scalar_matching::is_active()) return scalar_t(0.0);

    LOG_INFO("CQ2_susy::NLO 1");
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1});

    // double lu = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 7);
    // double ld = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 8);
    double z = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 1);
    double aY = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 11);
    // double yt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", 1);
    double Q_match = src.get_val(ParameterType::WILSON, "EW_SCALE", 1);;
    double mW = src.get_val(ParameterType::SM, "MASS", 24);
    double mH = src.get_val(ParameterType::BSM, "MASS", 37);
    double sw2 = src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4);
    double mass_b_muW = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1});
    double mass_top_muW = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 6);
    double g2 = src.get_val(ParameterType::SM, "GAUGE", 2);
    double tanb = src.get_val(ParameterType::BSM, "MINPAR", 3);

    double muQ = src.get_val(ParameterType::BSM, "HMIX", 1);
    // double xH = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 2);
    // double xH0 = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 3);
    // double alpha = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 9);
    // double beta = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 6);

    double kappa = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 19);
    double epsfac = src.get_val(ParameterType::WILSON, "EPSILON_SUSY", 5);

    complex_t NQ11H=-src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot)*(tanb)*tanb/4./mW/mW*(f141(xt,z)+8.*xt*(f30(xt,z)+xt*(f30(xt*1.0001,z)-f30(xt*0.9999,z))/0.0002)*log(Q_match*Q_match/mass_top_muW/mass_top_muW));
    complex_t BQ11H=src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot)*(tanb)*tanb/4./mW/mW*(f111(xt,z)+8.*(f70(xt*1.0001,z)-f70(xt*0.9999,z))/0.0002*log(Q_match*Q_match/mass_top_muW/mass_top_muW));
    complex_t CQ1H_1=(NQ11H+BQ11H)*mass_b_muW/sw2;
    complex_t CQ2H_1=-CQ1H_1;

    LOG_INFO("CQ2_susy::NLO 2");
    
    complex_t BQ11c1=0.;
    complex_t BQ11c2=0.;
    complex_t NQ11c=0.;
    complex_t NQ21c=0.;
    complex_t BQ11f1=0.;
    complex_t BQ11f2=0.;
    complex_t NQ11f=0.;
    complex_t NQ21f=0.;

    complex_t Dp{0}, Dm{0}, temp{0}, temp2{0};
    complex_t a0a{0}, a0b{0}, a0c{0}, a0Q1{0}, a0Q2{0}, a0p{0}, a1{0},a2p{0};

    for(int ie=0;ie<2;ie++) {
        for(int ae=0;ae<6;ae++){
            for(int je=0;je<2;je++)  {
                for(int be=0;be<3;be++) {
                    BQ11c1+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, je, ae, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6, ie, be, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5, je, be, 1})*(f121(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))+4.*(f50(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*1.0001,pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))-f50(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*0.9999,pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)))/0.0002*log(pow(mass_top_muW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}),2.))));
                    BQ11c2+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, je, ae, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5, ie, be, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6, je, be, 1})*fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}))*(f131(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))+4.*(f60(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*1.0001,pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))-f60(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*0.9999,pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)))/0.0002*log(pow(mass_top_muW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}),2.))));

                    for(int me=0;me<6;me++){ 
                        for(int ne=0;ne<3;ne++) {
                            Dp=0.;
                            Dm=0.;
                            for(int fe=1;fe<3;fe++) { 	
                                Dp+=src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe})/sqrt(2.)/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*muQ*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, me, fe})+src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, me, fe}));
                                Dm+=src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe})/sqrt(2.)/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*muQ*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, me, fe})-src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, me, fe}));
                            }
                            a0a=-(fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}))*(f181(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.))+4.*(f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.)*1.0001)-f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.)*0.9999))/0.0002*log(pow(mass_top_muW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}),2.)))*src.get_val(ParameterType::BSM, "UMIX", {ie+1, 1+1})*src.get_val(ParameterType::BSM, "VMIX", {je+1, 0+1}))*kron(ae,me);
                            a0b=-((f191(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.))+4.*(f40(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.)*1.0001)-f40(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.)*0.9999))/0.0002*log(pow(mass_top_muW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}),2.)))*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*src.get_val(ParameterType::BSM, "UMIX", {ie+1, 0+1}))*kron(ae,me);
                            a0c=1./mW*(f171(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, me})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))+4.*(f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, me})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))+(f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*1.0001,pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, me})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))-f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*0.9999,pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, me})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)))/0.0002+(f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, me})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*1.0001)-f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, me})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*0.9999))/0.0002)*log(pow(mass_top_muW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}),2.)))*kron(ie,je);
                        
                            a0Q1=a0a+a0b+Dp*a0c;
                            a0Q2=-a0a+a0b+Dm*a0c; //TODO : careful with RECKM -> VCKM (check with superiso values for real part)
                            a0p=4.*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne})/mW/(src.get_val(ParameterType::SM, "VCKM", {be, 2})*src.get_val(ParameterType::SM, "VCKM", {ne, 1})/src.get_val(ParameterType::SM, "VCKM", {2, 2})/src.get_val(ParameterType::SM, "VCKM", {2, 1}))/src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*f151(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))*kron(ie,je)*kron(ae,me)*kron(be,ne);
                            a1=src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/sqrt(2.)/mW*(f161(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))+4.*(f80(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*1.0001)-f80(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*0.9999))/0.0002*log(pow(mass_top_muW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}),2.)))*kron(ie,je)*kron(ae,me);
                            a2p=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, me, be})*(src.get_val(ParameterType::SM, "VCKM", {be, 2})*src.get_val(ParameterType::SM, "VCKM", {ne, 1})/src.get_val(ParameterType::SM, "VCKM", {2, 2})/src.get_val(ParameterType::SM, "VCKM", {2, 1}))*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})/2./mW*f151(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))*kron(ie,je)*kron(ae,me)*kron(be,ne);
                            
                            NQ11c+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, me, be})*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*(a0Q1+a1*tanb)
                            +src.get_val(ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne})*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*a0p
                            +src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, me, be})*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*a2p*pow(src.get_val(ParameterType::SM, "MASS", 3)*tanb,2.);	
                            NQ21c+=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, me, be})*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*(a0Q2+a1*tanb)
                            +src.get_val(ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne})*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*a0p
                            +src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, me, be})*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*a2p*pow(src.get_val(ParameterType::SM, "MASS", 3)*tanb,2.);
                        }
                    }
                
                }
                for(int be=0;be<6;be++) {
                    for(int ce=0;ce<6;ce++) {
                        for(int fe=0;fe<3;fe++) {
                            BQ11f1+=-src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, je, be, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2})*pow(mW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,ae, ce})*src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9, ce, be})*(1.+log(pow(mass_top_muW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}),2.)))	*(f90(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, fe})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6,ie, fe, 2})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5, je, fe, 2}));
                            BQ11f2+=-src.get_val(ParameterType::WILSON, "MATRIX_BSM", {3, je, be, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2})*pow(mW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.)*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,ae, ce})*src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9, ce, be})*(1.+log(pow(mass_top_muW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}),2.)))	*(fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}))*f100(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, be})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {16, fe})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.))*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {5, ie, fe, 1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {6, je, fe, 1}));

                        }
                    }
                }
            }

            LOG_INFO("CQ2_susy::NLO 3");

            for(int me=0;me<3;me++) {
                for(int ne=0;ne<3;ne++) {
                    for(int de=0;de<6;de++) {
                        for(int ke=0;ke<6;ke++) {
                            
                            temp2 = src.get_val(ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, me, ne})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, de, me})*src.get_val(ParameterType::BSM, "UMIX", {ie+1, 1+1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,ae, ke})*src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9, ke, de})*(1.+log(pow(mass_top_muW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}),2.)))*tanb*src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/sqrt(2.)*f30(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.));
                            NQ11f+=temp2;
                            NQ21f+=temp2;
                            for(int ce=0;ce<6;ce++) {
                                Dp=0.;
                                Dm=0.;
                                for(int fe=0;fe<3;fe++) 
                                {		
                                    Dp+=src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe})/sqrt(2.)/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*muQ*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, de, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ce, fe})+src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, de, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, ce, fe})); 
                                    Dm+=src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe})/sqrt(2.)/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*muQ*(src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, de, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ce, fe})-src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, de, fe})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {2, ce, fe})); 
                                }
                                temp=src.get_val(ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, me, ne})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, ce, me})*src.get_val(ParameterType::BSM, "UMIX", {ie+1, 1+1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,ae, ke})*src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9, ke, de})*
                                (1.+log(pow(mass_top_muW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}),2.)))*f60(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}),2.));	
                                NQ11f+=Dp*temp;
                                NQ21f+=Dm*temp;
                            }
                            for(int je=0;je<2;je++) {
                                temp=-src.get_val(ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, me, ne})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {1, de, me})*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9,ae, ke})*src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je})*src.get_val(ParameterType::WILSON, "MATRIX_BSM", {9, ke, de})*
                                (1.+log(pow(mass_top_muW/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}),2.)))*mW*(fabs(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}))*f60(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.))*src.get_val(ParameterType::BSM, "UMIX", {ie+1, 1+1})*src.get_val(ParameterType::BSM, "VMIX", {je+1, 0+1})+
                                f50(pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.),pow(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {14, de})/src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}),2.))*src.get_val(ParameterType::BSM, "UMIX", {je+1, 1+1})*src.get_val(ParameterType::BSM, "VMIX", {ie+1, 0+1})); 
                                NQ11f+=temp;
                                NQ21f+=-temp;
                            }

                        }
                    }
                }
                    
            }

            LOG_INFO("CQ2_susy::NLO 4");
        }
    }
    
    LOG_INFO("CQ2_susy::NLO 5");
    // complex_t BQ11c=(BQ11c1+BQ11c2)*kappa*mW*mW/2./g2/g2/sw2;
    complex_t BQ21c=-(BQ11c1-BQ11c2)*kappa*mW*mW/2./g2/g2/sw2;


    NQ11c*=src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot)*(tanb)*tanb/mW/(mH*mH-mW*mW)*aY*mass_b_muW/sw2;
    NQ21c*=-src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot)*(tanb)*tanb/mW/(mH*mH-mW*mW)*aY*mass_b_muW/sw2;
    
    complex_t CQ2charg_1=NQ21c+BQ21c;
        
    // complex_t BQ11f=(BQ11f1+BQ11f2)*2./3.*kappa/g2/g2/sw2;
    complex_t BQ21f=-(BQ11f1-BQ11f2)*2./3.*kappa/g2/g2/sw2;
    
    
    NQ11f*=-4./3.*src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot)*(tanb)*tanb/mW/mW/(mH*mH-mW*mW)*aY*mass_b_muW/sw2;

    NQ21f*=4./3.*src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot)*(tanb)*tanb/mW/mW/(mH*mH-mW*mW)*aY*mass_b_muW/sw2;

    complex_t CQ2four_1=NQ21f+BQ21f;

    complex_t coeff_temp = (CQ2H_1+CQ2charg_1)/epsfac+CQ2four_1;

    LOG_INFO("CQ2_susy::NLO 6");

    return coeff_temp;
}

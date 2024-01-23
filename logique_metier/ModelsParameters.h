#ifndef MODEL_PARAMETERS_H
#define MODEL_PARAMETERS_H

class Parameters {
public:
    double SM, gp, g2, MSOFT_Q, mass_top_pole, mass_b_pole;
    double A_b, tan_beta, mu_Q, mass_gluino, mass_b1, mass_b2;
    double inv_alpha_em, M2_Q, mass_t1, mass_t2, A_t, MqL3_Q, MbR_Q, mass_stl, mass_cha1, mass_cha2;
    std::vector<double> yut, yub, mass_neut;
    std::vector<std::vector<double>> sbot_mix, charg_Umix, charg_Vmix, stop_mix, neut_mix;
    std::vector<std::vector<double>> stop_tan_betamix;

    QCDParameters run;
    Parameters(); // Constructeur pour initialiser les paramètres
};

#endif // MODEL_PARAMETERS_H

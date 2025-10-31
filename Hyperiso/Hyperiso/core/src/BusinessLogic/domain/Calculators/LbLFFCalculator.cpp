#include "LbLFFCalculator.h"

LbLFFCalculator::LbLFFCalculator(LbL_FF_Src src) {
    ObsParameterProxy p(ParameterType::FLAVOR);
    double m_B = p("FMASS", 521);
    double m_K = p("FMASS", 311);
    this->m_Lb = p("FMASS", 5122);
    this->m_L = p("FMASS", 3122);
    this->t_p = std::pow(m_B + m_K, 2);
    this->t_0 = std::pow(this->m_Lb - this->m_L, 2);
    this->load_FF_params(src);
}

complex_t LbLFFCalculator::z(double t, double t_p, double t_0) {
    double a = std::sqrt(t_p - t);
    double b = std::sqrt(t_p - t_0);
    return (a - b) / (a + b);
}

double LbLFFCalculator::get(LbL_FF a, double q2) {
    return F_a(a, q2);
}

void LbLFFCalculator::load_FF_params(LbL_FF_Src src) {
    int ff_id = (int)(src) + 1;
    int sse_order = this->sse_order.at(src);
    std::string src_block = "Lb_L";
    
    auto get_m = [ff_id, src_block] (int i) { return ObsParameterProxy()(ParamId{ParameterType::DECAY, src_block, {ff_id, 0, i}}); };
    this->m_R[LbL_FF::F_PLUS] = this->m_R[LbL_FF::F_PERP] = this->m_R[LbL_FF::H_PLUS] = this->m_R[LbL_FF::H_PERP] = get_m(1);
    this->m_R[LbL_FF::G_PLUS] = this->m_R[LbL_FF::G_PERP] = this->m_R[LbL_FF::H_TILDE_PLUS] = this->m_R[LbL_FF::H_TILDE_PERP] = get_m(2);

    for (int i = 1; i <= 8; i++) {
        for (int j = 0; j <= sse_order; j++) {
            ParamId PId {ParameterType::DECAY, src_block, {ff_id, i, j}};
            this->alpha_ai[(LbL_FF)(i - 1)][j] = ObsParameterProxy()(PId);
        }
    }
}

double LbLFFCalculator::pole(double q2, double m_R) {
    return 1. / (1 - q2 / std::pow(m_R, 2));
}

double LbLFFCalculator::F_a(LbL_FF a, double q2) {
    auto ai = this->alpha_ai.at(a);
    double P = pole(q2, this->m_R.at(a));
    double Z = std::real(z(q2, this->t_p, this->t_0));
    return P * (ai[0] + Z * (ai[1] + Z * ai[2]));
}

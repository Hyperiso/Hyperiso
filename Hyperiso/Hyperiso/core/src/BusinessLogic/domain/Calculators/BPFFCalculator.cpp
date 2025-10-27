#include "BPFFCalculator.h"

BPFFCalculator::BPFFCalculator(int B_id, int P_id, BP_FF_Src src) {
    if (!this->allowed_decays.contains({B_id, P_id})) {
        LOG_ERROR("ValueError", "Wrong meson PDG code in BPFFCalculator constructor:", B_id, ",", P_id);
    }

    ObsParameterProxy p(ParameterType::FLAVOR);
    this->m_B = p("FMASS", B_id);
    this->m_P = p("FMASS", P_id);
    this->t_p = std::pow(this->m_B + this->m_P, 2);
    this->t_m = std::pow(this->m_B - this->m_P, 2);
    this->t_0 = src == BP_FF_Src::HPQCD22 ? 0. : this->t_p * (1. - std::sqrt(1 - this->t_m / this->t_p));

    if (src == BP_FF_Src::AS || src == BP_FF_Src::FLAG24 || src == BP_FF_Src::HPQCD22)
        this->z_0 = 0.0;
    else
        this->z_0 = std::real(z(0.0, this->t_p, this->t_0));

    this->src_block = allowed_decays.at({B_id, P_id});
    this->load_FF_params(src);
}

complex_t BPFFCalculator::z(double t, double t_p, double t_0) {
    double a = std::sqrt(t_p - t);
    double b = std::sqrt(t_p - t_0);
    return (a - b) / (a + b);
}

double BPFFCalculator::get(BP_FF a, double q2) {
    switch (a) {
    case BP_FF::XI_P:
        return this->F_a(BP_FF::F_PLUS, q2);
    default:
        return this->F_a(a, q2);
    }
}

void BPFFCalculator::load_FF_params(BP_FF_Src src) {
    int ff_id = (int)(src) + 1;
    int order = this->sse_order.at(src);
    std::string src_block = this->src_block;
    
    auto get_m = [ff_id, src_block] (int i) { return ObsParameterProxy()(ParamId{ParameterType::DECAY, src_block, {ff_id, 0, i}}); };
    this->m_R[BP_FF::F_PLUS] = get_m(1);
    this->m_R[BP_FF::F_0] = get_m(2);
    this->m_R[BP_FF::F_T] = get_m(3);

    for (int i = 1; i <= 3; i++) {
        for (int j = 0; j <= order; j++) {
            ParamId PId {ParameterType::DECAY, src_block, {ff_id, i, j}};
            this->alpha_ai[(BP_FF)(i - 1)][j] = ObsParameterProxy()(PId);
        }
    }

    this->L_chi = 1.0;

    if (src == BP_FF_Src::AS) {
        this->alpha_ai[BP_FF::F_0][3] = (2 * this->alpha_ai[BP_FF::F_0][2] - this->alpha_ai[BP_FF::F_0][1]) / 3;
        this->alpha_ai[BP_FF::F_T][3] = (2 * this->alpha_ai[BP_FF::F_T][2] - this->alpha_ai[BP_FF::F_T][1]) / 3;
        this->m_R[BP_FF::F_PLUS] += m_B;
        this->m_R[BP_FF::F_T] += m_B;
    }

    if (src == BP_FF_Src::FLAG24) {
        this->alpha_ai[BP_FF::F_PLUS][3] = (2 * this->alpha_ai[BP_FF::F_PLUS][2] - this->alpha_ai[BP_FF::F_PLUS][1]) / 3;
        double z0 = std::real(this->z(0, this->t_p, this->t_0));
        this->alpha_ai[BP_FF::F_0][2] = this->alpha_ai[BP_FF::F_PLUS][2] + (this->alpha_ai[BP_FF::F_PLUS][0] - this->alpha_ai[BP_FF::F_0][0]) / pow(z0, 2.) + (this->alpha_ai[BP_FF::F_PLUS][1] - this->alpha_ai[BP_FF::F_0][1]) / z0 + this->alpha_ai[BP_FF::F_PLUS][3] * z0;
    }

    if (src == BP_FF_Src::HPQCD22) {
        this->alpha_ai[BP_FF::F_PLUS][3] = (2 * this->alpha_ai[BP_FF::F_PLUS][2] - this->alpha_ai[BP_FF::F_PLUS][1]) / 3;
        this->alpha_ai[BP_FF::F_T][3] = (2 * this->alpha_ai[BP_FF::F_T][2] + this->alpha_ai[BP_FF::F_T][1]) / 3;
        this->L_chi = ObsParameterProxy(ParameterType::DECAY)(this->src_block, 16);
    }
}

double BPFFCalculator::pole(double q2, double m_R) {
    if (fpeq(m_R, 0.0)) return 1.0;

    return 1. / (1 - q2 / std::pow(m_R, 2));
}

double BPFFCalculator::E(double q2) {
    return (this->m_B * this->m_B + this->m_P * this->m_P - q2) / (2 * this->m_B);
}

double BPFFCalculator::F_a(BP_FF a, double q2) {
    auto ai = this->alpha_ai.at(a);
    double P = pole(q2, this->m_R.at(a));
    double Z = std::real(z(q2, this->t_p, this->t_0)) - this->z_0;

    double res {0.0};
    for (int k = 0; k < 4; k++) {
        res += ai[k] * std::pow(Z, k);
    }

    return this->L_chi * P * res;
}
#include "BVFFCalculator.h"

BVFFCalculator::BVFFCalculator(int B_id, int V_id, BV_FF_Src src) {
    if (!this->allowed_decays.contains({B_id, V_id})) {
        LOG_ERROR("ValueError", "Wrong meson PDG code in BVFFCalculator constructor:", B_id, ",", V_id);
    }

    if (V_id == 333 && (src == BV_FF_Src::GKvD_SR || src == BV_FF_Src::GKvD_SR_LAT)) {
        LOG_WARN("GKvD formfactors are not available for Bs > phi decays. Defaulting to BFS formfactors.");
        src = src == BV_FF_Src::GKvD_SR ? BV_FF_Src::BSZ_SR : BV_FF_Src::BSZ_SR_LAT;
    }

    ObsParameterProxy p(ParameterType::FLAVOR);
    this->m_B = p("FMASS", B_id);
    this->m_B2 = std::pow(this->m_B, 2);
    this->m_B4 = std::pow(this->m_B2, 2);
    this->m_V = p("FMASS", V_id);
    this->m_V2 = std::pow(this->m_V, 2);
    this->m_V4 = std::pow(this->m_V2, 2);
    this->t_p = std::pow(this->m_B + this->m_V, 2);
    this->t_m = std::pow(this->m_B - this->m_V, 2);
    this->t_0 = src == BV_FF_Src::HLMW ? 12. : this->t_p * (1. - std::sqrt(1 - this->t_m / this->t_p));
    this->z_0 = std::real(z(0.0, this->t_p, this->t_0));
    this->src_block = allowed_decays.at({B_id, V_id});
    this->load_FF_params(src);
}

complex_t BVFFCalculator::z(double t, double t_p, double t_0) {
    double a = std::sqrt(t_p - t);
    double b = std::sqrt(t_p - t_0);
    return (a - b) / (a + b);
}

double BVFFCalculator::get(BV_FF a, double q2) {
    switch (a) {
    case BV_FF::A2:
        return this->A_2(q2);
    case BV_FF::T3:
        return this->T_3(q2);
    case BV_FF::XI_PERP:
        return this->xi_perp(q2);
    case BV_FF::XI_PAR:
        return this->xi_par(q2);
    case BV_FF::F_PERP:
        return this->f_perp(q2);
    case BV_FF::F_PAR:
        return this->f_par(q2);
    case BV_FF::F_0:
        return this->f_0(q2);
    default:
        return this->F_a(a, q2);
    }
}

void BVFFCalculator::load_FF_params(BV_FF_Src src) {
    int ff_id = (int)(src) + 1;
    int sse_order = src == BV_FF_Src::HLMW ? 1 : 2;
    std::string src_block = this->src_block;
    
    auto get_m = [ff_id, src_block] (int i) { return ObsParameterProxy()(ParamId{ParameterType::DECAY, src_block, {ff_id, 0, i}}); };
    this->m_R[BV_FF::A0] = get_m(1);
    this->m_R[BV_FF::V] = this->m_R[BV_FF::T1] = get_m(2);
    this->m_R[BV_FF::A1] = this->m_R[BV_FF::A12] = this->m_R[BV_FF::T2] = this->m_R[BV_FF::T23] = get_m(3);

    for (int i = 1; i <= 7; i++) {
        for (int j = 0; j <= sse_order; j++) {
            ParamId PId {ParameterType::DECAY, src_block, {ff_id, i, j}};
            this->alpha_ai[(BV_FF)(i - 1)][j] = ObsParameterProxy()(PId);
        }
    }
}

double BVFFCalculator::pole(double q2, double m_R) {
    return 1. / (1 - q2 / std::pow(m_R, 2));
}

double BVFFCalculator::E(double q2) {
    return (this->m_B2 + this->m_V2 - q2) / (2 * this->m_B);
}

double BVFFCalculator::F_a(BV_FF a, double q2) {
    auto ai = this->alpha_ai.at(a);
    double P = pole(q2, this->m_R.at(a));
    double Z = std::real(z(q2, this->t_p, this->t_0)) - this->z_0;
    return P * (ai[0] + Z * (ai[1] + Z * ai[2]));
}

double BVFFCalculator::A_2(double q2) {
    double A_1 = F_a(BV_FF::A1, q2);
    double A_12 = F_a(BV_FF::A12, q2);
    return (this->t_p * (this->m_B2 - this->m_V2 - q2) * A_1 - 16. * this->m_B * this->m_V2 * (this->m_B + this->m_V) * A_12) / ((this->t_p - q2) * (this->t_m - q2));
}

double BVFFCalculator::T_3(double q2) {
    double T_2 = F_a(BV_FF::T2, q2);
    double T_23 = F_a(BV_FF::T23, q2);
    return ((this->m_B2 - this->m_V2) * (this->m_B2 + 3. * this->m_V2 - q2) * T_2 - 8. * this->m_B * this->m_V2 * (this->m_B - this->m_V) * T_23) / ((this->t_p - q2) * (this->t_m - q2));
}

double BVFFCalculator::xi_perp(double q2) {
    return this->m_B * F_a(BV_FF::V, q2) / (this->m_B + this->m_V);
}

double BVFFCalculator::xi_par(double q2) {
    return (this->m_B + this->m_V) * F_a(BV_FF::A1, q2) / (2. * E(q2)) - (this->m_B - this->m_V) * A_2(q2) / this->m_B;
}

double BVFFCalculator::f_perp(double q2) {
    double lambda = this->m_B4 + this->m_V4 + q2 * q2 - 2 * (this->m_B2 * this->m_V2 + q2 * (this->m_B2 + this->m_V2));
    return std::sqrt(2. * lambda) / (this->m_B + this->m_V) * F_a(BV_FF::V, q2);
}

double BVFFCalculator::f_par(double q2) {
    return RT2 * (this->m_B + this->m_V) * F_a(BV_FF::A1, q2);
}

double BVFFCalculator::f_0(double q2) {
    double lambda = this->m_B4 + this->m_V4 + q2 * q2 - 2 * (this->m_B2 * this->m_V2 + q2 * (this->m_B2 + this->m_V2));
    return ((this->m_B * this->m_B - q2 - this->m_V * this->m_V) * this->t_p * F_a(BV_FF::A1, q2) - lambda * A_2(q2)) / (2. * this->m_V * (this->m_B + this->m_V) * sqrt(q2));
}
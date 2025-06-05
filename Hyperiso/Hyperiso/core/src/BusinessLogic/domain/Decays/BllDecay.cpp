#include "BllDecay.h"


scalar_t BllDecay::W1(scalar_t r, scalar_t CQ1, scalar_t CPQ1, bool prime) {
    CPQ1 = prime ? CPQ1 : scalar_t(0.);
    return r * (CQ1 - CPQ1);
}

scalar_t BllDecay::W2Q(scalar_t r, scalar_t CQ2, scalar_t CPQ2, bool prime) {
    CPQ2 = prime ? CPQ2 : scalar_t(0.);
    return r * (CQ2 - CPQ2);
}

scalar_t BllDecay::W210(scalar_t x, scalar_t C10, scalar_t CP10, bool prime) {
    CP10 = prime ? CP10 : scalar_t(0.);
    return 2. * (C10 - CP10) * x;
}

scalar_t BllDecay::ckm(scalar_t V_tb, scalar_t V_tq) {
    return std::pow(std::abs(V_tb * std::conj(V_tq)), 2);
}

scalar_t BllDecay::BR_avg_Bq_mumu(scalar_t w1,
                                scalar_t w2q,
                                scalar_t w210,
                                scalar_t ckm,
                                scalar_t x,
                                scalar_t G_F,
                                scalar_t inv_alpha,
                                scalar_t f_B,
                                scalar_t m_B,
                                scalar_t life_B)
{
    scalar_t b = std::sqrt(1. - 4. * std::pow(x, 2));
    scalar_t pref = std::pow(G_F * f_B / inv_alpha, 2) * std::pow(m_B / M_PI, 3) * life_B * ckm * b / (64. * HBAR);
    std::cout << "b : " << b << std::endl;
    std::cout << "pref : " << pref << std::endl;
    std::cout << "w1: " << w1 << std::endl;
    std::cout << "w2q: " << w2q << std::endl;
    std::cout << "w210: " << w210 << std::endl;

    return pref * (std::pow(b * std::abs(w1), 2) + std::pow(std::abs(w2q + w210), 2));
}

scalar_t BllDecay::A_DG(scalar_t x, scalar_t r, scalar_t w210, scalar_t w1q, scalar_t w2q, scalar_t C10_SM) {
    std::cout << C10_SM << std::endl;
    scalar_t S = x * std::sqrt(1. - 4. * x * x) * w1q / (2. * C10_SM);
    scalar_t P = (w210 / x + w2q * x) / (2. * C10_SM);
    scalar_t abs_S = std::pow(std::abs(S), 2);
    scalar_t abs_P = std::pow(std::abs(P), 2);

    return (abs_P * std::cos(2 * std::arg(P)) - abs_S * std::cos(2 * std::arg(S))) / (abs_P + abs_S);
}

scalar_t BllDecay::BR_untag_Bs_mumu(scalar_t br_avg, scalar_t ys, scalar_t A) {
    scalar_t untag_factor = (1. + A * ys) / (1. - ys * ys);
    return untag_factor * br_avg;
}

void BllDecay::build_op_tree() {
    
    // SM Parameters
    auto inv_alpha_em = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "SMINPUTS", 1));
    auto G_F = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "SMINPUTS", 2));

    auto m_mu = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 13));
    auto m_d = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 1));
    auto m_s = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 3));
    auto V_tb = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "VCKM", LhaID(2, 2)));
    auto V_ts = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "VCKM", LhaID(2, 1)));
    auto V_td = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "VCKM", LhaID(2, 0)));
    auto m_b_pole = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "QCD", LhaID(5, 2)));

    // Wilson Coefficients
    auto C10 = std::make_shared<ParameterNode>(ParamId(ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), LhaID(3051313, 4137, 2)));
    auto C10_SM = std::make_shared<ParameterNode>(ParamId(ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), LhaID(3051313, 4137, 0)));
    auto CQ1 = std::make_shared<ParameterNode>(ParamId(ParameterType::WILSON, GroupMapper::str(WGroup::BScalar, ScaleType::HADRONIC), LhaID(3051313, 3230, 2)));
    auto CQ2 = std::make_shared<ParameterNode>(ParamId(ParameterType::WILSON, GroupMapper::str(WGroup::BScalar, ScaleType::HADRONIC), LhaID(3051313, 3233, 2)));
    auto CP10 = std::make_shared<ParameterNode>(ParamId(ParameterType::WILSON, GroupMapper::str(WGroup::BPrime, ScaleType::HADRONIC), LhaID(3051313, 4234, 2)));
    auto CPQ1 = std::make_shared<ParameterNode>(ParamId(ParameterType::WILSON, GroupMapper::str(WGroup::BPrime, ScaleType::HADRONIC), LhaID(3051313, 3130, 2)));
    auto CPQ2 = std::make_shared<ParameterNode>(ParamId(ParameterType::WILSON, GroupMapper::str(WGroup::BPrime, ScaleType::HADRONIC), LhaID(3051313, 3133, 2)));

    // Flavor Parameters
    auto m_Bs = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 531));
    auto life_Bs = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FLIFE", 531));
    auto f_Bs = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FCONST", LhaID(531,1)));
    auto m_Bd = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 511));
    auto life_Bd = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FLIFE", 511));
    auto f_Bd = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FCONST", LhaID(511,1)));

    // Misc experimental input
    auto y_s = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_ll", 1)); // y_s = life_Bs * Delta(Gamma_s) / 2
    
    // Operator nodes
    auto xs = std::make_shared<OperatorNode>("xs", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] / values[1]; });
    xs->addChildren({m_mu, m_Bs});
    auto rs = std::make_shared<OperatorNode>("rs", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] / (values[1] + values[2]); });
    rs->addChildren({m_Bs, m_b_pole, m_s});
    auto w1s = std::make_shared<OperatorNode>("W1s", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->W1(values[0], values[1], values[2], true); });
    w1s->addChildren({rs, CQ1, CPQ1});
    auto w2qs = std::make_shared<OperatorNode>("W2Qs", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->W2Q(values[0], values[1], values[2], true); });
    w2qs->addChildren({rs, CQ2, CPQ2});
    auto w210s = std::make_shared<OperatorNode>("W210s", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->W210(values[0], values[1], values[2], true); });
    w210s->addChildren({xs, C10, CP10});
    auto ckm_s = std::make_shared<OperatorNode>("CKM_s", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->ckm(values[0], values[1]); });
    ckm_s->addChildren({V_tb, V_ts});
    auto br_avg_Bs_mumu = std::make_shared<OperatorNode>("BR_Bs__mu_mu", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->BR_avg_Bq_mumu(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9]); });
    br_avg_Bs_mumu->addChildren({w1s, w2qs, w210s, ckm_s, xs, G_F, inv_alpha_em, f_Bs, m_Bs, life_Bs});
    roots.emplace(Observables::BR_BS_MUMU, br_avg_Bs_mumu);

    auto xd = std::make_shared<OperatorNode>("xd", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] / values[1]; });
    xd->addChildren({m_mu, m_Bd});
    auto rd = std::make_shared<OperatorNode>("rd", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] / (values[1] + values[2]); });
    rd->addChildren({m_Bd, m_b_pole, m_d});
    //TODO or not TODO : why was it addChild with only rd ? not addChildren with CQ1, CPQ1, CQ2, etc. ?
    auto w1d = std::make_shared<OperatorNode>("W1d", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->W1(values[0], values[1], values[2], false); });
    w1d->addChildren({rd, CQ1, CPQ1});
    auto w2qd = std::make_shared<OperatorNode>("W2Qd", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->W2Q(values[0], values[1], values[2], false); });
    w2qd->addChildren({rd, CQ2, CPQ2});
    auto w210d = std::make_shared<OperatorNode>("W210d", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->W210(values[0], values[1], values[2], false); });
    w210d->addChildren({xd, C10, CP10});
    auto ckm_d = std::make_shared<OperatorNode>("CKM_d", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->ckm(values[0], values[1]); });
    ckm_d->addChildren({V_tb, V_td});

    auto br_avg_Bd_mumu = std::make_shared<OperatorNode>("BR_Bd__mu_mu", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->BR_avg_Bq_mumu(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9]); });
    br_avg_Bd_mumu->addChildren({w1d, w2qd, w210d, ckm_d, xd, G_F, inv_alpha_em, f_Bd, m_Bd, life_Bd});
    roots.emplace(Observables::BR_BD_MUMU, br_avg_Bd_mumu);
    auto a_dg = std::make_shared<OperatorNode>("A_DeltaGamma", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->A_DG(values[0], values[1], values[2], values[3], values[4], values[5]); });
    a_dg->addChildren({xs, rs, w210s, w1s, w2qs, C10_SM});
    auto br_untag_Bs_mumu = std::make_shared<OperatorNode>("BR_untag_Bs__mu_mu", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->BR_untag_Bs_mumu(values[0], values[1], values[2]); });
    br_untag_Bs_mumu->addChildren({br_avg_Bs_mumu, y_s, a_dg});
    roots.emplace(Observables::BR_BS_MUMU_UNTAG, br_untag_Bs_mumu);
}

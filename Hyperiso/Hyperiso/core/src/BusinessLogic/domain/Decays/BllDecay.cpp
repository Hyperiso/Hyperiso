#include "BllDecay.h"


complex_t BllDecay::W1(double r, bool prime) {
    complex_t cq1 = w_proxy->getFR(WGroup::BScalar, WCoef::CQ1, QCDOrder::NLO);
    complex_t cpq1 = prime ? w_proxy->getFR(WGroup::BPrime, WCoef::CPQ1, QCDOrder::NLO) : 0;
    return r * (cq1 - cpq1);
}

complex_t BllDecay::W2Q(double r, bool prime) {
    complex_t cq2 = w_proxy->getFR(WGroup::BScalar, WCoef::CQ2, QCDOrder::NLO);
    complex_t cpq2 = prime ? w_proxy->getFR(WGroup::BPrime, WCoef::CPQ2, QCDOrder::NLO) : 0;
    return r * (cq2 - cpq2);
}

complex_t BllDecay::W210(double x, bool prime) {
    complex_t c10 = w_proxy->getFR(WGroup::B, WCoef::C10, QCDOrder::NNLO);
    complex_t cp10 = prime ? w_proxy->getFR(WGroup::BPrime, WCoef::CP10, QCDOrder::NLO) : 0;
    return 2. * (c10 - cp10) * x;
}

double BllDecay::ckm(complex_t V_tb, complex_t V_tq) {
    return std::pow(std::abs(V_tb * std::conj(V_tq)), 2);
}

double BllDecay::BR_avg_Bq_mumu(complex_t w1,
                                complex_t w2q,
                                complex_t w210,
                                double ckm,
                                double x,
                                double G_F,
                                double inv_alpha,
                                double f_B,
                                double m_B,
                                double life_B)
{
    double b = std::sqrt(1 - 4 * std::pow(x, 2));
    double pref = std::pow(G_F * f_B / inv_alpha, 2) * std::pow(m_B / M_PI, 3) * life_B * ckm * b / (64 * HBAR);
    return pref * (std::pow(b * std::abs(w1), 2) + std::pow(std::abs(w2q + w210), 2));
}

double BllDecay::A_DG(double x, double r) {
    complex_t cq1 = w_proxy->getFR(WGroup::BScalar, WCoef::CQ1, QCDOrder::NLO);
    complex_t cpq1 = w_proxy->getFR(WGroup::BPrime, WCoef::CPQ1, QCDOrder::LO);
    complex_t cq2 = w_proxy->getFR(WGroup::BScalar, WCoef::CQ2, QCDOrder::NLO);
    complex_t cpq2 = w_proxy->getFR(WGroup::BPrime, WCoef::CPQ2, QCDOrder::LO);
    complex_t c10 = w_proxy->getFR(WGroup::B, WCoef::C10, QCDOrder::NNLO);
    complex_t cp10 = w_proxy->getFR(WGroup::BPrime, WCoef::CP10, QCDOrder::LO);
    complex_t c10sm = w_proxy->getFR(WGroup::B, WCoef::C10, QCDOrder::NNLO, true);

    complex_t S = x * std::sqrt(1 - 4 * x * x) * r * (cq1 - cpq1) / (2. * c10sm);
    complex_t P = (c10 - cp10) / c10sm + x * r * (cq2 - cpq2) / (2. * c10sm);
    double abs_S = std::pow(std::abs(S), 2);
    double abs_P = std::pow(std::abs(P), 2);

    return (abs_P * std::cos(2 * std::arg(P)) - abs_S * std::cos(2 * std::arg(S))) / (abs_P + abs_S);
}

double BllDecay::BR_untag_Bs_mumu(double br_avg, double ys, double A) {
    double untag_factor = (1 + A * ys) / (1 - ys * ys);
    return untag_factor * br_avg;
}

void BllDecay::build_op_tree() {
    
    // SM Parameters
    auto inv_alpha_em = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "SMINPUTS", 1));
    auto G_F = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "SMINPUTS", 2));

    auto m_mu = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 13));
    auto m_d = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 1));
    auto m_s = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 3));
    auto V_tb = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "RECKM", 22));
    auto V_ts = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "RECKM", 21));
    auto V_td = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "RECKM", 20));
    auto m_b_pole = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "QCD", LhaID(5, 2)));

    // Flavor Parameters
    auto m_Bs = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 531));
    auto life_Bs = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FLIFE", 531));
    auto f_Bs = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FCONST", 53101));
    auto m_Bd = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 511));
    auto life_Bd = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FLIFE", 511));
    auto f_Bd = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FCONST", 51101));

    // Misc experimental input
    auto y_s = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_ll", 1)); // y_s = life_Bs * Delta(Gamma_s) / 2
    
    // Operator nodes
    auto xs = std::make_shared<OperatorNode>("xs", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] / values[1]; });
    xs->addChildren({m_mu, m_Bs});
    auto rs = std::make_shared<OperatorNode>("rs", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] / (values[1] + values[2]); });
    rs->addChildren({m_Bs, m_b_pole, m_s});
    auto w1s = std::make_shared<OperatorNode>("W1s", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->W1(values[0], true); });
    w1s->addChildren({rs});
    auto w2qs = std::make_shared<OperatorNode>("W2Qs", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->W2Q(values[0], true); });
    w2qs->addChildren({rs});
    auto w210s = std::make_shared<OperatorNode>("W210s", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->W210(values[0], true); });
    w210s->addChildren({xs});
    auto ckm_s = std::make_shared<OperatorNode>("CKM_s", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->ckm(values[0], values[1]); });
    ckm_s->addChildren({V_tb, V_ts});
    auto br_avg_Bs_mumu = std::make_shared<OperatorNode>("BR_Bs__mu_mu", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->BR_avg_Bq_mumu(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9]); });
    br_avg_Bs_mumu->addChildren({w1s, w2qs, w210s, ckm_s, xs, G_F, inv_alpha_em, f_Bs, m_Bs, life_Bs});
    roots.emplace(Observables::BR_BS_MUMU, br_avg_Bs_mumu);

    auto xd = std::make_shared<OperatorNode>("xd", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] / values[1]; });
    xd->addChildren({m_mu, m_Bd});
    auto rd = std::make_shared<OperatorNode>("rd", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] / (values[1] + values[2]); });
    rd->addChildren({m_Bd, m_b_pole, m_d});
    auto w1d = std::make_shared<OperatorNode>("W1d", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->W1(values[0], false); });
    w1d->addChild(rd);
    auto w2qd = std::make_shared<OperatorNode>("W2Qd", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->W2Q(values[0], false); });
    w2qd->addChild(rd);
    auto w210d = std::make_shared<OperatorNode>("W210d", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->W210(values[0], false); });
    w210d->addChild(xd);
    auto ckm_d = std::make_shared<OperatorNode>("CKM_d", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->ckm(values[0], values[1]); });
    ckm_d->addChildren({V_tb, V_td});
    auto br_avg_Bd_mumu = std::make_shared<OperatorNode>("BR_Bd__mu_mu", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->BR_avg_Bq_mumu(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9]); });
    br_avg_Bd_mumu->addChildren({w1d, w2qd, w210d, ckm_d, xd, G_F, inv_alpha_em, f_Bd, m_Bd, life_Bd});
    roots.emplace(Observables::BR_BD_MUMU, br_avg_Bd_mumu);

    auto a_dg = std::make_shared<OperatorNode>("A_DeltaGamma", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->A_DG(values[0], values[1]); });
    a_dg->addChildren({xs, rs});
    auto br_untag_Bs_mumu = std::make_shared<OperatorNode>("BR_untag_Bs__mu_mu", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->BR_untag_Bs_mumu(values[0], values[1], values[2]); });
    br_untag_Bs_mumu->addChildren({br_avg_Bs_mumu, y_s, a_dg});
    roots.emplace(Observables::BR_BS_MUMU_UNTAG, br_untag_Bs_mumu);
}

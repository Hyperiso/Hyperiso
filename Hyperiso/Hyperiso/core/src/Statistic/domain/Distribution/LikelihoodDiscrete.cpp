#include "LikelihoodDiscrete.h"

LikelihoodDiscrete::LikelihoodDiscrete(std::vector<double> values,
                                       std::vector<double> weights,
                                       unsigned int seed,
                                       bool standardize)
    : eng_(seed), u01_(0.0, 1.0), values_(std::move(values)), standardize_(standardize)
{
    if (values_.empty()) throw std::invalid_argument("LikelihoodDiscrete: empty values.");
    if (weights.size() != values_.size()) throw std::invalid_argument("LikelihoodDiscrete: weights size mismatch.");

    // compute mean/std if standardize requested (using normalized weights)
    double wsum = 0.0;
    for (double w : weights) {
        if (!(w >= 0.0) || !std::isfinite(w)) throw std::invalid_argument("LikelihoodDiscrete: invalid weight (must be finite >=0).");
        wsum += w;
    }
    if (wsum <= 0.0) throw std::invalid_argument("LikelihoodDiscrete: sum(weights) must be > 0.");

    if (standardize_) {
        // normalize weights for moments
        double m = 0.0;
        for (std::size_t i = 0; i < values_.size(); ++i) m += values_[i] * (weights[i] / wsum);

        double v = 0.0;
        for (std::size_t i = 0; i < values_.size(); ++i) {
            double d = values_[i] - m;
            v += d * d * (weights[i] / wsum);
        }
        mean_ = m;
        std_  = (v > 0.0) ? std::sqrt(v) : 1.0;
    }

    build_alias_tables(std::move(weights));
}

void LikelihoodDiscrete::build_alias_tables(std::vector<double> weights) {
    const std::size_t n = values_.size();
    prob_.assign(n, 0.0);
    alias_.assign(n, 0);

    // Normalize weights to average 1 (scaled by n)
    double sumw = std::accumulate(weights.begin(), weights.end(), 0.0);
    if (sumw <= 0.0) throw std::invalid_argument("LikelihoodDiscrete: sum(weights) must be > 0.");

    std::vector<double> scaled(n);
    for (std::size_t i = 0; i < n; ++i) scaled[i] = (weights[i] * n) / sumw;

    std::vector<std::size_t> small;
    std::vector<std::size_t> large;
    small.reserve(n);
    large.reserve(n);

    for (std::size_t i = 0; i < n; ++i) {
        if (scaled[i] < 1.0) small.push_back(i);
        else large.push_back(i);
    }

    while (!small.empty() && !large.empty()) {
        const std::size_t s = small.back(); small.pop_back();
        const std::size_t l = large.back(); large.pop_back();

        prob_[s] = scaled[s];
        alias_[s] = l;

        scaled[l] = (scaled[l] + scaled[s]) - 1.0;
        if (scaled[l] < 1.0) small.push_back(l);
        else large.push_back(l);
    }

    for (std::size_t i : large) {
        prob_[i] = 1.0;
        alias_[i] = i;
    }
    for (std::size_t i : small) {
        prob_[i] = 1.0;
        alias_[i] = i;
    }
}

Vector LikelihoodDiscrete::sample(std::size_t n) {
    Vector out(n);
    const std::size_t m = values_.size();
    std::uniform_int_distribution<std::size_t> uid(0, m - 1);

    for (std::size_t k = 0; k < n; ++k) {
        const std::size_t i = uid(eng_);
        const double r = u01_(eng_);
        const std::size_t idx = (r < prob_[i]) ? i : alias_[i];
        double x = values_[idx];
        if (standardize_) x = (x - mean_) / std_;
        out[k] = x;
    }
    return out;
}

#include "JointDistribution.h"

JointDistribution::JointDistribution(
    std::vector<std::unique_ptr<IMarginalDistribution>> marginals,
    std::unique_ptr<ICopula> copula) :
    marginals_(std::move(marginals)), copula_(std::move(copula))
{}

std::vector<Vector> JointDistribution::sample(std::size_t n) const {
    std::vector<Vector> u = this->copula_->sample_u(n);
    std::vector<Vector> x (n, Vector(u[0].size(), 0.0));

    for (size_t i = 0; i < marginals_.size(); i++) {
        for (size_t j = 0; j < n; j++) {
            x[j][i] = marginals_.at(i)->ppf(u[j][i]);
        }
    }
    
    return x;
}

Vector JointDistribution::sample() const {
    return sample(1)[0];
}

double JointDistribution::logpdf(Vector x) const {
    if (x.size() != marginals_.size()) 
        throw std::invalid_argument("Wrong size of random vector.");

    Vector u = Vector(x.size(), 0.0);
    double log_marg {0.0};

    for (size_t i = 0; i < marginals_.size(); i++) {
        u[i] = marginals_[i]->cdf(x[i]);
        log_marg += marginals_[i]->logpdf(x[i]);
    }

    return log_marg + copula_->log_density(u);    
}

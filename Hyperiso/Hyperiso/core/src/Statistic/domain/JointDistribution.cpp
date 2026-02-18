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
    Vector u = this->copula_->sample_u();
    Vector x (u.size(), 0.0);

    for (size_t i = 0; i < marginals_.size(); i++) {
        x[i] = marginals_.at(i)->ppf(u[i]);
    }
    
    return x;
}

double JointDistribution::logpdf(Vector x) const {
    if (x.size() != marginals_.size()) 
        throw std::invalid_argument("Wrong size of random vector.");

    Vector u = Vector(x.size(), 0.0);
    double log_marg {0.0};

    // printf("size(x) = %i\n", x.size());

    for (size_t i = 0; i < marginals_.size(); i++) {
        u[i] = std::clamp(marginals_[i]->cdf(x[i]), 1e-13, 1 - 1e-13);
        log_marg += marginals_[i]->logpdf(x[i]);
        // printf("x_i = %.5e\n", x[i]);
        // printf("marginals_i mean = %.5e\n", marginals_[i]->mean());
        // printf("marginals_i std = %.5e\n", marginals_[i]->std());
        // printf("log_marg_i = %.5e\n", marginals_[i]->logpdf(x[i]));
    }

    // printVector(x);
    // printf("log_marg = %.5e\n", log_marg);
    // printf("log_copula = %.5e\n", copula_->log_density(u));

    return log_marg + copula_->log_density(u);    
}

std::size_t JointDistribution::dim() {
    return marginals_.size();
}

Vector JointDistribution::get_stds() {
    Vector stds;
    for (auto& m : this->marginals_)
        stds.emplace_back(m->std());

    return stds;
}

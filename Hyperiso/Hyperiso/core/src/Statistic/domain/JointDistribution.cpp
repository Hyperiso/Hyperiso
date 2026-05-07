#include "JointDistribution.h"

JointDistribution::JointDistribution(
    std::vector<std::unique_ptr<IMarginalDistribution>> marginals,
    std::unique_ptr<ICopula> copula) :
    marginals_(std::move(marginals)), copula_(std::move(copula))
{}

std::vector<std::vector<double>> JointDistribution::sample(std::size_t n) const {
    std::vector<std::vector<double>> u = this->copula_->sample_u(n);
    std::vector<std::vector<double>> x (n, std::vector<double>(u[0].size(), 0.0));

    for (size_t i = 0; i < marginals_.size(); i++) {
        for (size_t j = 0; j < n; j++) {
            x[j][i] = marginals_.at(i)->ppf(u[j][i]);
        }
    }
    
    return x;
}

std::vector<double> JointDistribution::sample() const {
    std::vector<double> u = this->copula_->sample_u();
    std::vector<double> x (u.size(), 0.0);

    for (size_t i = 0; i < marginals_.size(); i++) {
        x[i] = marginals_.at(i)->ppf(u[i]);
    }
    
    return x;
}

double JointDistribution::logpdf(std::vector<double> x) const {
    if (x.size() != marginals_.size()) 
        throw std::invalid_argument("Wrong size of random vector.");

    std::vector<double> u = std::vector<double>(x.size(), 0.0);
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

    // printstd::vector<double>(x);
    // printf("log_marg = %.5e\n", log_marg);
    // printf("log_copula = %.5e\n", copula_->log_density(u));

    return log_marg + copula_->log_density(u);    
}

RealMatrix JointDistribution::curvature(std::vector<double> x) const {
    if (x.size() != marginals_.size()) 
        throw std::invalid_argument("Wrong size of random vector.");

    std::size_t d = x.size();

    std::vector<double> u = std::vector<double>(d, 0.0);
    RealMatrix W(d, d);
    double marginal_term {0.0};
    Vector f_i (d, 0.0);
    Vector df_i (d, 0.0);
    Vector ddf_i (d, 0.0);

    for (size_t i = 0; i < marginals_.size(); i++) {
        u[i] = std::clamp(marginals_[i]->cdf(x[i]), 1e-13, 1 - 1e-13);
        PDFDiff fdf = marginals_[i]->f_df_ddf(x[i]);
        f_i[i] = fdf.f;
        df_i[i] = fdf.df;
        ddf_i[i] = fdf.ddf;
    }

    LogDensityDiff cdc = copula_->log_c_dc_ddc(u);

    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            W.at(i, j) = -cdc.ddlog_c.at(i, j) * f_i[i] * f_i[j];
            if (i == j) W.at(i, j) += -cdc.dlog_c.at(i, 0) * df_i[i] - (ddf_i[i] / f_i[i] - std::pow(df_i[i] / f_i[i], 2));
        }
    }

    return W;  
}

std::size_t JointDistribution::dim() {
    return marginals_.size();
}

std::vector<double> JointDistribution::get_stds() {
    std::vector<double> stds;
    for (auto& m : this->marginals_)
        stds.emplace_back(m->std());

    return stds;
}

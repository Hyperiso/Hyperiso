#include "WithGaussianConstraints.h"

WithGaussianConstraints::WithGaussianConstraints(
    std::shared_ptr<ILikelihood> base,
    std::shared_ptr<JointDistribution> constraints_dist,
    std::vector<std::size_t> constrained_params)
    :   base(std::move(base)), 
        constraints_dist(std::move(constraints_dist)), 
        constrained_params(std::move(constrained_params))
{}

double WithGaussianConstraints::nll(const Vector &theta) const {
    double base_nll = this->base->nll(theta);

    std::vector<double> constrained_vals;
    for (std::size_t i : this->constrained_params) {
        constrained_vals.emplace_back(theta[i]);
    }
    
    double constraint_nll = -this->constraints_dist->logpdf(constrained_vals);

    return base_nll + constraint_nll;
}

std::vector<fit_app::ParameterDefinition> WithGaussianConstraints::get_param_defs() const {
    return base->get_param_defs();
}

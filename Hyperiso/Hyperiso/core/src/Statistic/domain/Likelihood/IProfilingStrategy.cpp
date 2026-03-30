#include "IProfilingStrategy.h"

IProfilingStrategy::IProfilingStrategy(std::size_t x_id, std::size_t y_id, const FitResult& fr)
    : x_id(x_id), y_id(y_id), fr(fr)
{}

ProfileRequest SliceProfilingStrategy::build_request(
    double px,
    double py,
    const std::map<std::size_t, double>& current_argmin
) const {
    const std::size_t dim_p = fr.p_hat.size();
    const std::size_t dim_nuis = fr.eta_hat.size();
    const std::size_t dim = dim_p + dim_nuis;

    ProfileRequest pr;
    pr.start.resize(dim);

    for (std::size_t i = 0; i < dim_p; ++i) {
        double val = fr.p_hat[i];
        if (i == x_id) val = px;
        if (i == y_id) val = py;
        pr.fixed_params[i] = val;
        pr.start[i] = val;
    }

    for (std::size_t i = 0; i < dim_nuis; ++i) {
        const std::size_t idx = dim_p + i;
        pr.free_params.push_back(idx);

        auto it = current_argmin.find(idx);
        if (it == current_argmin.end()) {
            throw std::runtime_error("SliceProfilingStrategy: missing warm-start entry for nuisance index "
                                     + std::to_string(idx));
        }
        pr.start[idx] = it->second;
    }

    return pr;
}


std::map<std::size_t, double> SliceProfilingStrategy::init_warm_start() const {
    std::map<std::size_t, double> start;
    const std::size_t dim_p = fr.p_hat.size();
    const std::size_t dim_nuis = fr.eta_hat.size();

    for (std::size_t i = 0; i < dim_nuis; ++i) {
        start[dim_p + i] = fr.eta_hat.at(i);
    }

    return start;
}

ProfileRequest ProjectionProfilingStrategy::build_request(
    double px,
    double py,
    const std::map<std::size_t, double>& current_argmin
) const {
    const std::size_t dim_p = fr.p_hat.size();
    const std::size_t dim_nuis = fr.eta_hat.size();
    const std::size_t dim = dim_p + dim_nuis;

    ProfileRequest pr;
    pr.start.resize(dim);

    for (std::size_t i = 0; i < dim; ++i) {
        if (i == x_id) {
            pr.fixed_params[i] = px;
            pr.start[i] = px;
        } else if (i == y_id) {
            pr.fixed_params[i] = py;
            pr.start[i] = py;
        } else {
            pr.free_params.push_back(i);

            auto it = current_argmin.find(i);
            if (it == current_argmin.end()) {
                throw std::runtime_error("ProjectionProfilingStrategy: missing warm-start entry for free index "
                                         + std::to_string(i));
            }
            pr.start[i] = it->second;
        }
    }

    return pr;
}

std::map<std::size_t, double> ProjectionProfilingStrategy::init_warm_start() const {
    std::map<std::size_t, double> start;
    const std::size_t dim_p = fr.p_hat.size();
    const std::size_t dim_nuis = fr.eta_hat.size();

    for (std::size_t i = 0; i < dim_p + dim_nuis; ++i) {
        start[i] = (i < dim_p) ? fr.p_hat.at(i) : fr.eta_hat.at(i - dim_p);
    }

    start.erase(x_id);
    start.erase(y_id);

    return start;
}
#include "IProfilingStrategy.h"

IProfilingStrategy::IProfilingStrategy(std::size_t x_id, std::size_t y_id, const FitResult& fr)
    : x_id(x_id), y_id(y_id), fr(fr)
{}

ProfileRequest SliceProfilingStrategy::build_request(double px, double py) const
{
    std::size_t dim_p = fr.p_hat.size();
    std::size_t dim_nuis = fr.eta_hat.size();
    std::size_t dim = dim_p + dim_nuis;
    std::map<std::size_t, double> start = this->get_start();

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
        std::size_t idx = dim_p + i;
        pr.free_params.push_back(idx);
        pr.start[idx] = start.at(i);
    }

    return pr;
}


std::map<std::size_t, double> SliceProfilingStrategy::get_start() const {
    std::map<std::size_t, double> start;
    std::size_t dim_p = fr.p_hat.size();
    std::size_t dim_nuis = fr.eta_hat.size();

    for (size_t i = 0; i < dim_nuis; i++) {
        start[dim_p + i] = fr.eta_hat.at(i);
    }

    return start;
}

ProfileRequest ProjectionProfilingStrategy::build_request(double px, double py) const
{
    std::size_t dim_p = fr.p_hat.size();
    std::size_t dim_nuis = fr.eta_hat.size();
    std::size_t dim = dim_p + dim_nuis;
    std::map<std::size_t, double> start = this->get_start();

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
            pr.start[i] = start.at(i);
        }
    }

    return pr;
}

std::map<std::size_t, double> ProjectionProfilingStrategy::get_start() const {
    std::map<std::size_t, double> start;
    std::size_t dim_p = fr.p_hat.size();
    std::size_t dim_nuis = fr.eta_hat.size();

    for (size_t i = 0; i < dim_p + dim_nuis; i++) {
        start[i] = i < dim_p ? fr.p_hat.at(i) : fr.eta_hat.at(i - dim_p);
    }

    start.erase(x_id);
    start.erase(y_id);

    return start;
}

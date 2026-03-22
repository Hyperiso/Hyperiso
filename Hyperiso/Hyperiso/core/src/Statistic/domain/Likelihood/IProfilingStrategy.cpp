#include "IProfilingStrategy.h"

IProfilingStrategy::IProfilingStrategy(std::size_t x_id, std::size_t y_id, const FitResult& fr)
    : x_id(x_id), y_id(y_id), fr(fr)
{}

// ProfileRequest SliceProfilingStrategy::build_request(
//     double px, double py,
//     const std::map<std::size_t, double> &current_argmin) const
// {
//     std::size_t dim_p = fr.p_hat.size();
//     std::size_t dim_nuis = fr.eta_hat.size();
//     std::size_t dim = dim_p + dim_nuis;

//     ProfileRequest pr;

//     for (size_t i = 0; i < dim; i++) {
//         if (i < dim_p)
//             pr.fixed_params.emplace(i, fr.p_hat[i]);
//         else
//             pr.free_params.emplace_back(i);
//     }

//     pr.fixed_params[x_id] = px;
//     pr.fixed_params[y_id] = py;

//     auto p_start = unzip(pr.fixed_params).vals;
//     auto nuis_start = unzip(current_argmin).vals;
//     p_start.insert(p_start.end(), nuis_start.begin(), nuis_start.end());
//     pr.start = p_start;

//     return pr;
// }
//TODO : Niels : change that because start(dim) + emplace back -> bad size
// pr.start = unzip(current_argmin).vals loose fixed value
ProfileRequest SliceProfilingStrategy::build_request(
    double px, double py,
    const std::map<std::size_t, double>& current_argmin) const
{
    std::size_t dim_p = fr.p_hat.size();
    std::size_t dim_nuis = fr.eta_hat.size();
    std::size_t dim = dim_p + dim_nuis;

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
        pr.start[idx] = current_argmin.at(idx);
    }

    return pr;
}


std::map<std::size_t, double> SliceProfilingStrategy::init_warm_start() const {
    std::map<std::size_t, double> start;
    std::size_t dim_p = fr.p_hat.size();
    std::size_t dim_nuis = fr.eta_hat.size();

    for (size_t i = 0; i < dim_nuis; i++) {
        start[dim_p + i] = fr.eta_hat.at(i);
    }

    return start;
}

// ProfileRequest ProjectionProfilingStrategy::build_request(
//     double px, double py,
//     const std::map<std::size_t, double> &current_argmin) const
// {
//     std::size_t dim_p = fr.p_hat.size();
//     std::size_t dim_nuis = fr.eta_hat.size();
//     std::size_t dim = dim_p + dim_nuis;
//     std::vector<double> start (dim);

//     ProfileRequest pr;

//     for (size_t i = 0; i < dim; i++) {
//         if (i == x_id)
//             start.emplace_back(px);
//         else if (i == y_id)
//             start.emplace_back(py);
//         else {
//             start.emplace_back(current_argmin.at(i));
//             pr.free_params.emplace_back(i);
//         }
//     }

//     pr.fixed_params[x_id] = px;
//     pr.fixed_params[y_id] = py;

//     pr.start = unzip(current_argmin).vals;

//     return pr;
// }

ProfileRequest ProjectionProfilingStrategy::build_request(
    double px, double py,
    const std::map<std::size_t, double>& current_argmin) const
{
    std::size_t dim_p = fr.p_hat.size();
    std::size_t dim_nuis = fr.eta_hat.size();
    std::size_t dim = dim_p + dim_nuis;

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
            pr.start[i] = current_argmin.at(i);
        }
    }

    return pr;
}

std::map<std::size_t, double> ProjectionProfilingStrategy::init_warm_start() const {
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

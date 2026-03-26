#ifndef __IPROFILINGSTRATEGY_H__
#define __IPROFILINGSTRATEGY_H__

#include "Profiler.h"

struct FitResult {
    std::vector<double> p_hat; // MLE estimators
    std::vector<double> eta_hat; // profiled-at-MLE nuisances
    std::vector<double> p_hat_std;
    RealMatrix p_hat_correlations;
    double ell_hat {0.0}; // minimum NLL
};

class IProfilingStrategy {
public:
    IProfilingStrategy(std::size_t x_id, std::size_t y_id, const FitResult& fr);

    virtual ~IProfilingStrategy() = default;
    virtual ProfileRequest build_request(double px, double py) const = 0;
    virtual std::map<std::size_t, double> get_start() const = 0;

    std::size_t get_x_id() { return x_id; };
    std::size_t get_y_id() { return y_id; };

protected:
    std::size_t x_id, y_id;
    FitResult fr;
};

class SliceProfilingStrategy : public IProfilingStrategy {
public:
    SliceProfilingStrategy(std::size_t x_id, std::size_t y_id, const FitResult& fr) : IProfilingStrategy(x_id, y_id, fr) {}
    ProfileRequest build_request(double px, double py) const override;
    std::map<std::size_t, double> get_start() const override;
};

class ProjectionProfilingStrategy : public IProfilingStrategy {
public:
    ProjectionProfilingStrategy(std::size_t x_id, std::size_t y_id, const FitResult& fr) : IProfilingStrategy(x_id, y_id, fr) {}
    ProfileRequest build_request(double px, double py) const override;
    std::map<std::size_t, double> get_start() const override;
};

#endif // __IPROFILINGSTRATEGY_H__

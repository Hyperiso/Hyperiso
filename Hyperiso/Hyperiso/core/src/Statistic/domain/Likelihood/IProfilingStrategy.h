#ifndef __IPROFILINGSTRATEGY_H__
#define __IPROFILINGSTRATEGY_H__

#include "Profiler.h"
#include "Fit.h"

class IProfilingStrategy {
public:
    IProfilingStrategy(std::size_t x_id, std::size_t y_id, const FitResult& fr);

    virtual ~IProfilingStrategy() = default;
    virtual ProfileRequest build_request(double px, double py, const std::map<std::size_t, double>& current_argmin) const = 0;

protected:
    std::size_t x_id, y_id;
    FitResult fr;
};

class SliceProfilingStrategy : public IProfilingStrategy {
public:
    SliceProfilingStrategy(std::size_t x_id, std::size_t y_id, const FitResult& fr) : IProfilingStrategy(x_id, y_id, fr) {}
    ProfileRequest build_request(double px, double py, const std::map<std::size_t, double>& current_argmin) const override;
};

class ProjectionProfilingStrategy : public IProfilingStrategy {
public:
    ProjectionProfilingStrategy(std::size_t x_id, std::size_t y_id, const FitResult& fr) : IProfilingStrategy(x_id, y_id, fr) {}
    ProfileRequest build_request(double px, double py, const std::map<std::size_t, double>& current_argmin) const override;
};

#endif // __IPROFILINGSTRATEGY_H__

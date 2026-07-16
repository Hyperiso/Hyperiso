#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Matrix.h"

namespace fit_app {

struct ParameterDefinition {
    std::string name;
    double value = 0.0;
    double step_hint = 0.0;
    std::optional<std::pair<double, double>> limits;
    bool fixed = false;
};

struct FitOptions {
    double up = 0.5;
    unsigned strategy = 2;
    unsigned max_fcn = 100000;
    double tolerance = 0.2;
    bool run_hesse = true;
    unsigned hesse_maxcalls = 0;
    bool verbose = true;
};

struct ContourOptionsBackEnd {
    double up = 0.5;
    unsigned npoints = 80;
    unsigned strategy = 2;
    unsigned max_fcn = 30000;
    double tolerance = 0.2;
};

// struct FitDiagnostics {
//     double fmin = std::numeric_limits<double>::quiet_NaN();
//     double edm = std::numeric_limits<double>::quiet_NaN();
//     int nfcn = -1;

//     bool ok = false;
//     bool has_valid_covar = false;
//     bool has_posdef_covar = false;
//     bool has_accurate_covar = false;
//     bool made_posdef = false;

//     std::vector<double> cov_eigs;
//     double cond_number = std::numeric_limits<double>::infinity();
// };

struct FitDiagnostics {
    double fmin = std::numeric_limits<double>::quiet_NaN();
    double edm = std::numeric_limits<double>::quiet_NaN();
    int nfcn = -1;

    bool ok = false;

    bool has_valid_parameters = false;
    bool has_valid_covar = false;
    bool has_posdef_covar = false;
    bool has_accurate_covar = false;
    bool made_posdef = false;
    bool hesse_failed = false;

    bool reached_call_limit = false;
    bool above_max_edm = false;

    std::vector<double> cov_eigs;
    double cond_number = std::numeric_limits<double>::infinity();
};

class BackendState {
public:
    virtual ~BackendState() = default;
};

struct BackendFitResult {
    std::vector<double>values;
    std::vector<double> errors;
    RealMatrix covariance;
    FitDiagnostics diagnostics;
    std::shared_ptr<const BackendState> state;
};

struct BackendContourResult {
    bool success = false;
    std::vector<std::pair<double, double>> points;
};

class IObjectiveFunction {
public:
    virtual ~IObjectiveFunction() = default;
    virtual double operator()(const std::vector<double>& x) const = 0;
    virtual double error_definition() const = 0;
};

class LambdaObjectiveFunction final : public IObjectiveFunction {
public:
    LambdaObjectiveFunction(std::function<double(const std::vector<double>&)> fn, double up)
        : fn_(std::move(fn)), up_(up) {}

    double operator()(const std::vector<double>& x) const override {
        return fn_(x);
    }

    double error_definition() const override {
        return up_;
    }

private:
    std::function<double(const std::vector<double>&)> fn_;
    double up_ = 0.5;
};

class IFitBackend {
public:
    virtual ~IFitBackend() = default;

    virtual BackendFitResult minimize(const IObjectiveFunction& objective,
                                      const std::vector<ParameterDefinition>& parameters,
                                      const FitOptions& options) const = 0;

    virtual BackendFitResult minimize_with_fixed(const IObjectiveFunction& objective,
                                                 const std::vector<ParameterDefinition>& parameters,
                                                 const FitOptions& options,
                                                 const std::vector<std::size_t>& fixed_indices,
                                                 const std::vector<double>& fixed_values) const = 0;

    virtual BackendContourResult contour(const IObjectiveFunction& objective,
                                         const BackendFitResult& reference_fit,
                                         std::size_t x_index,
                                         std::size_t y_index,
                                         const ContourOptionsBackEnd& options) const = 0;
};

std::unique_ptr<IFitBackend> make_minuit_backend();

template <class T>
std::string to_string_any(const T& x) {
    std::ostringstream oss;
    oss << x;
    return oss.str();
}

inline std::vector<double> linspace(double a, double b, std::size_t n) {
    std::vector<double> out(n);
    if (n == 0) return out;
    if (n == 1) {
        out[0] = a;
        return out;
    }
    for (std::size_t i = 0; i < n; ++i) {
        out[i] = a + (b - a) * double(i) / double(n - 1);
    }
    return out;
}

inline double safe_step(double value, double scale_hint) {
    const double abs_value = std::fabs(value);
    const double abs_scale = std::fabs(scale_hint);

    double step = 0.0;
    if (std::isfinite(abs_scale) && abs_scale > 0.0) step = 0.05 * abs_scale;
    if (std::isfinite(abs_value) && abs_value > 0.0) step = std::max(step, 0.01 * abs_value);
    if (!std::isfinite(step) || step <= 0.0) step = 1e-3;
    return step;
}

} // namespace fit_app

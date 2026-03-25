#include "FitAbstraction.h"

#include <algorithm>
#include <iostream>

#include "minuit-cpp/FCNBase.hh"
#include "minuit-cpp/FunctionMinimum.hh"
#include "minuit-cpp/MnContours.hh"
#include "minuit-cpp/MnEigen.hh"
#include "minuit-cpp/MnHesse.hh"
#include "minuit-cpp/MnMigrad.hh"
#include "minuit-cpp/MnUserParameters.hh"

namespace M2 = MinuitCpp;

namespace fit_app {
namespace {

class MinuitObjective final : public M2::FCNBase {
public:
    explicit MinuitObjective(const IObjectiveFunction& objective)
        : objective_(objective) {}

    double operator()(const std::vector<double>& x) const override {
        try {
            const double value = objective_(x);
            return std::isfinite(value) ? value : 1e300;
        } catch (...) {
            return 1e300;
        }
    }

    double Up() const override {
        return objective_.error_definition();
    }

private:
    const IObjectiveFunction& objective_;
};

class MinuitState final : public BackendState {
public:
    explicit MinuitState(M2::FunctionMinimum minimum)
        : minimum_(std::move(minimum)) {}

    const M2::FunctionMinimum& minimum() const {
        return minimum_;
    }

private:
    M2::FunctionMinimum minimum_;
};

class MinuitCppBackend final : public IFitBackend {
public:
    BackendFitResult minimize(const IObjectiveFunction& objective,
                              const std::vector<ParameterDefinition>& parameters,
                              const FitOptions& options) const override {
        return do_minimize(objective, parameters, options, {}, {});
    }

    BackendFitResult minimize_with_fixed(const IObjectiveFunction& objective,
                                         const std::vector<ParameterDefinition>& parameters,
                                         const FitOptions& options,
                                         const std::vector<std::size_t>& fixed_indices,
                                         const std::vector<double>& fixed_values) const override {
        if (fixed_indices.size() != fixed_values.size()) {
            throw std::invalid_argument("fixed_indices/fixed_values size mismatch");
        }
        return do_minimize(objective, parameters, options, fixed_indices, fixed_values);
    }

    BackendContourResult contour(const IObjectiveFunction& objective,
                                 const BackendFitResult& reference_fit,
                                 std::size_t x_index,
                                 std::size_t y_index,
                                 const ContourOptions& options) const override {
        BackendContourResult out;

        auto state = std::dynamic_pointer_cast<const MinuitState>(reference_fit.state);
        if (!state) return out;

        try {
            MinuitObjective minuit_objective(objective);
            M2::MnContours mn_contours(minuit_objective, state->minimum(), options.strategy);
            out.points = mn_contours(static_cast<unsigned>(x_index), static_cast<unsigned>(y_index), options.npoints);
            out.success = out.points.size() >= 4;
        } catch (...) {
            out.points.clear();
            out.success = false;
        }

        return out;
    }

private:
    static void apply_parameters(M2::MnUserParameters& upar,
                                 const std::vector<ParameterDefinition>& parameters) {
        for (const auto& parameter : parameters) {
            upar.Add(parameter.name.c_str(), parameter.value, safe_step(parameter.value, parameter.step_hint));
            if (parameter.limits.has_value()) {
                upar.SetLimits(parameter.name.c_str(), parameter.limits->first, parameter.limits->second);
            }
            if (parameter.fixed) {
                upar.Fix(parameter.name.c_str());
            }
        }
    }

    static void apply_fixed_values(M2::MnMigrad& migrad,
                                   const std::vector<std::size_t>& fixed_indices,
                                   const std::vector<double>& fixed_values) {
        for (std::size_t i = 0; i < fixed_indices.size(); ++i) {
            const unsigned idx = static_cast<unsigned>(fixed_indices[i]);
            migrad.SetValue(idx, fixed_values[i]);
            migrad.Fix(idx);
        }
    }

    static void log_summary(const M2::FunctionMinimum& minimum) {
        std::cout << "\n=== [MINUIT] FunctionMinimum ===\n";
        std::cout << "IsValid            = " << minimum.IsValid() << "\n";
        std::cout << "HasValidParameters = " << minimum.HasValidParameters() << "\n";
        std::cout << "HasValidCovariance = " << minimum.HasValidCovariance() << "\n";
        std::cout << "HasAccurateCovar   = " << minimum.HasAccurateCovar() << "\n";
        std::cout << "HasPosDefCovar     = " << minimum.HasPosDefCovar() << "\n";
        std::cout << "HasMadePosDefCovar = " << minimum.HasMadePosDefCovar() << "\n";
        std::cout << "HesseFailed        = " << minimum.HesseFailed() << "\n";
        std::cout << "Fval               = " << minimum.Fval() << "\n";
        std::cout << "EDM                = " << minimum.Edm() << "\n";
        std::cout << "NFcn               = " << minimum.NFcn() << "\n";
        std::cout << "Up                 = " << minimum.Up() << "\n";
    }

    static BackendFitResult extract_result(const M2::FunctionMinimum& minimum,
                                           const std::vector<ParameterDefinition>& parameters,
                                           bool extract_full_covariance) {
        BackendFitResult out;
        out.values.resize(parameters.size());
        out.errors.assign(parameters.size(), 0.0);
        out.covariance = RealMatrix(parameters.size(), parameters.size());
        out.state = std::make_shared<MinuitState>(minimum);

        out.diagnostics.fmin = minimum.Fval();
        out.diagnostics.edm = minimum.Edm();
        out.diagnostics.nfcn = minimum.NFcn();
        out.diagnostics.ok = minimum.IsValid();
        out.diagnostics.has_valid_parameters = minimum.HasValidParameters();
        out.diagnostics.hesse_failed = minimum.HesseFailed();
        out.diagnostics.has_valid_covar = minimum.HasValidCovariance();
        out.diagnostics.has_posdef_covar = minimum.HasPosDefCovar();
        out.diagnostics.has_accurate_covar = minimum.HasAccurateCovar();
        out.diagnostics.made_posdef = minimum.HasMadePosDefCovar();
        out.diagnostics.reached_call_limit = minimum.HasReachedCallLimit();
        out.diagnostics.above_max_edm = minimum.IsAboveMaxEdm();

        const auto& state = minimum.UserState();
        for (std::size_t i = 0; i < parameters.size(); ++i) {
            out.values[i] = state.Value(parameters[i].name.c_str());
            out.errors[i] = parameters[i].fixed ? 0.0 : state.Error(parameters[i].name.c_str());
        }

        // if (!extract_full_covariance || !minimum.HasValidCovariance()) {
        //     out.diagnostics.has_valid_covar = false;
        //     out.diagnostics.has_posdef_covar = false;
        //     out.diagnostics.has_accurate_covar = false;
        //     out.diagnostics.made_posdef = false;
        //     return out;
        // }
        if (!minimum.HasValidCovariance()) {
            return out;
        }

        if (!extract_full_covariance) {
            return out;
        }
        const auto& cov = state.Covariance();
        for (std::size_t i = 0; i < parameters.size(); ++i) {
            for (std::size_t j = 0; j < parameters.size(); ++j) {
                out.covariance.at(i, j) = cov(i, j);
            }
        }

        M2::MnEigen eigen;
        out.diagnostics.cov_eigs = eigen(cov);

        double min_pos = std::numeric_limits<double>::infinity();
        double max_pos = 0.0;
        for (double eig : out.diagnostics.cov_eigs) {
            if (std::isfinite(eig) && eig > 0.0) {
                min_pos = std::min(min_pos, eig);
                max_pos = std::max(max_pos, eig);
            }
        }

        if (min_pos < std::numeric_limits<double>::infinity() && max_pos > 0.0) {
            out.diagnostics.cond_number = max_pos / min_pos;
        }

        return out;
    }

    static BackendFitResult do_minimize(const IObjectiveFunction& objective,
                                        const std::vector<ParameterDefinition>& parameters,
                                        const FitOptions& options,
                                        const std::vector<std::size_t>& fixed_indices,
                                        const std::vector<double>& fixed_values) {
        MinuitObjective minuit_objective(objective);

        M2::MnUserParameters upar;
        apply_parameters(upar, parameters);

        M2::MnMigrad migrad(minuit_objective, upar, options.strategy);
        apply_fixed_values(migrad, fixed_indices, fixed_values);

        M2::FunctionMinimum minimum = migrad(options.max_fcn, options.tolerance);

        if (options.run_hesse) {
            M2::MnHesse hesse(options.strategy);
            hesse(minuit_objective, minimum, options.hesse_maxcalls);
        }

        if (options.verbose) {
            log_summary(minimum);
        }

        const bool has_explicitly_fixed_parameters = !fixed_indices.empty() ||
            std::any_of(parameters.begin(), parameters.end(), [](const ParameterDefinition& parameter) {
                return parameter.fixed;
            });

        return extract_result(minimum, parameters, !has_explicitly_fixed_parameters);
    }
};

} // namespace

std::unique_ptr<IFitBackend> make_minuit_backend() {
    return std::make_unique<MinuitCppBackend>();
}

} // namespace fit_app

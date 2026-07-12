#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>

#include <array>
#include <memory>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "StatisticInterface.h"

#include "AbstractConfig.h"
#include "MarginalType.h"
#include "IMarginalDistribution.h"
#include "GaussianMarginal.h"
#include "SplitGaussianMarginal.h"
#include "FlatMarginal.h"
#include "LikelihoodMarginal.h"
#include "MarginalFactory.h"
#include "MarginalConfigFactory.h"

#include "CopulaType.h"
#include "ICopula.h"
#include "GaussianCopula.h"
#include "StudentTCopula.h"
#include "CopulaFactory.h"
#include "JointDistribution.h"

#include "ILikelihood.h"
#include "BaseLikelihood.h"
#include "Fit.h"
#include "MCEngine.h"
#include "IProfilingStrategy.h"

#if __has_include("Profiler.h")
#  include "Profiler.h"
#  define STAT_BIND_HAS_PROFILER 1
#elif __has_include("WithProfiling.h")
#  include "WithProfiling.h"
#  define STAT_BIND_HAS_PROFILER 1
#else
#  define STAT_BIND_HAS_PROFILER 0
#endif

#if __has_include("WithGaussianConstraints.h")
#  include "WithGaussianConstraints.h"
#  define STAT_BIND_HAS_GAUSSIAN_CONSTRAINTS 1
#else
#  define STAT_BIND_HAS_GAUSSIAN_CONSTRAINTS 0
#endif

namespace py = pybind11;

namespace {

unsigned int fresh_seed() {
    return std::random_device{}();
}

// void init_abstract_config(py::module_& m) {
//     // Utile si ce fichier est chargé seul. Si AbstractConfig est déjà exposée
//     // ailleurs dans le même module, on évite une double déclaration.
//     if (!py::hasattr(m, "AbstractConfig")) {
//         py::class_<AbstractConfig>(m, "AbstractConfig");
//     }
// }

MarginalConfig to_marginal_config(const py::handle& obj) {
    if (py::isinstance<FlatMarginalCfg>(obj)) {
        return obj.cast<FlatMarginalCfg>();
    }
    if (py::isinstance<GaussianMarginalCfg>(obj)) {
        return obj.cast<GaussianMarginalCfg>();
    }
    if (py::isinstance<SplitGaussianMarginalCfg>(obj)) {
        return obj.cast<SplitGaussianMarginalCfg>();
    }
    if (py::isinstance<LikelihoodMarginalCfg>(obj)) {
        return obj.cast<LikelihoodMarginalCfg>();
    }
    throw py::type_error(
        "Marginal config must be one of: FlatMarginalCfg, GaussianMarginalCfg, "
        "SplitGaussianMarginalCfg, LikelihoodMarginalCfg."
    );
}

CopulaConfig to_copula_config(const py::handle& obj) {
    if (py::isinstance<GaussianCopulaConfig>(obj)) {
        return obj.cast<GaussianCopulaConfig>();
    }
    if (py::isinstance<StudentTCopulaConfig>(obj)) {
        return obj.cast<StudentTCopulaConfig>();
    }
    throw py::type_error("Copula config must be GaussianCopulaConfig or StudentTCopulaConfig.");
}

void init_marginals(py::module_& m) {
    py::enum_<MarginalType>(m, "MarginalType")
        .value("GAUSSIAN", MarginalType::GAUSSIAN)
        .value("HALF_GAUSSIAN", MarginalType::HALF_GAUSSIAN)
        .value("FLAT", MarginalType::FLAT)
        .value("LIKELIHOOD", MarginalType::LIKELIHOOD)
        .export_values();

    py::class_<FlatMarginalCfg, AbstractConfig>(m, "FlatMarginalCfg")
        .def(py::init<double, double>(), py::arg("a"), py::arg("b"))
        .def_readwrite("a", &FlatMarginalCfg::a)
        .def_readwrite("b", &FlatMarginalCfg::b)
        .def("__repr__", [](const FlatMarginalCfg& c) {
            std::ostringstream oss;
            oss << "FlatMarginalCfg(a=" << c.a << ", b=" << c.b << ")";
            return oss.str();
        });

    py::class_<GaussianMarginalCfg, AbstractConfig>(m, "GaussianMarginalCfg")
        .def(py::init<>())
        .def(py::init<double, double>(), py::arg("mu"), py::arg("sigma"))
        .def_readwrite("mu", &GaussianMarginalCfg::mu)
        .def_readwrite("sigma", &GaussianMarginalCfg::sigma)
        .def("__repr__", [](const GaussianMarginalCfg& c) {
            std::ostringstream oss;
            oss << "GaussianMarginalCfg(mu=" << c.mu << ", sigma=" << c.sigma << ")";
            return oss.str();
        });

    py::class_<SplitGaussianMarginalCfg, AbstractConfig>(m, "SplitGaussianMarginalCfg")
        .def(py::init<>())
        .def(py::init([](double mu, double sigma_p, double sigma_m) {
            SplitGaussianMarginalCfg c;
            c.mu = mu;
            c.sigma_p = sigma_p;
            c.sigma_m = sigma_m;
            return c;
        }), py::arg("mu"), py::arg("sigma_p"), py::arg("sigma_m"))
        .def_readwrite("mu", &SplitGaussianMarginalCfg::mu)
        .def_readwrite("sigma_p", &SplitGaussianMarginalCfg::sigma_p)
        .def_readwrite("sigma_m", &SplitGaussianMarginalCfg::sigma_m)
        .def("__repr__", [](const SplitGaussianMarginalCfg& c) {
            std::ostringstream oss;
            oss << "SplitGaussianMarginalCfg(mu=" << c.mu
                << ", sigma_p=" << c.sigma_p
                << ", sigma_m=" << c.sigma_m << ")";
            return oss.str();
        });

    py::class_<LikelihoodMarginalCfg, AbstractConfig>(m, "LikelihoodMarginalCfg")
        .def(py::init<>())
        .def(py::init([](Vector values, Vector weights) {
            LikelihoodMarginalCfg c;
            c.values = std::move(values);
            c.weights = std::move(weights);
            return c;
        }), py::arg("values"), py::arg("weights"))
        .def_readwrite("values", &LikelihoodMarginalCfg::values)
        .def_readwrite("weights", &LikelihoodMarginalCfg::weights)
        .def("__repr__", [](const LikelihoodMarginalCfg& c) {
            std::ostringstream oss;
            oss << "LikelihoodMarginalCfg(values=" << c.values.size()
                << ", weights=" << c.weights.size() << ")";
            return oss.str();
        });

    py::class_<IMarginalDistribution, std::shared_ptr<IMarginalDistribution>>(m, "IMarginalDistribution")
        .def("rvs", &IMarginalDistribution::rvs, py::arg("n"))
        .def("logpdf", &IMarginalDistribution::logpdf, py::arg("x"))
        .def("cdf", &IMarginalDistribution::cdf, py::arg("x"))
        .def("ppf", &IMarginalDistribution::ppf, py::arg("p"))
        .def("mean", &IMarginalDistribution::mean)
        .def("std", &IMarginalDistribution::std);

    py::class_<GaussianMarginal, IMarginalDistribution, std::shared_ptr<GaussianMarginal>>(m, "GaussianMarginal")
        .def(py::init([](double mu, double sigma) {
            return std::make_shared<GaussianMarginal>(mu, sigma, fresh_seed());
        }), py::arg("mu"), py::arg("sigma"))
        .def(py::init<double, double, unsigned int>(),
             py::arg("mu"), py::arg("sigma"), py::arg("seed"));

    py::class_<SplitGaussianMarginal, IMarginalDistribution, std::shared_ptr<SplitGaussianMarginal>>(m, "SplitGaussianMarginal")
        .def(py::init([](double mu, double sigma_p, double sigma_m) {
            return std::make_shared<SplitGaussianMarginal>(mu, sigma_p, sigma_m, fresh_seed());
        }), py::arg("mu"), py::arg("sigma_p"), py::arg("sigma_m"))
        .def(py::init<double, double, double, unsigned int>(),
             py::arg("mu"), py::arg("sigma_p"), py::arg("sigma_m"), py::arg("seed"));

    py::class_<FlatMarginal, IMarginalDistribution, std::shared_ptr<FlatMarginal>>(m, "FlatMarginal")
        .def(py::init([](double a, double b) {
            return std::make_shared<FlatMarginal>(a, b, fresh_seed());
        }), py::arg("a"), py::arg("b"))
        .def(py::init<double, double, unsigned int>(),
             py::arg("a"), py::arg("b"), py::arg("seed"));

    py::class_<LikelihoodMarginal, IMarginalDistribution, std::shared_ptr<LikelihoodMarginal>>(m, "LikelihoodMarginal")
        .def(py::init([](Vector values, Vector weights, bool standardize) {
            return std::make_shared<LikelihoodMarginal>(
                std::move(values), std::move(weights), fresh_seed(), standardize
            );
        }), py::arg("values"), py::arg("weights"), py::arg("standardize") = false)
        .def(py::init<Vector, Vector, unsigned int, bool>(),
             py::arg("values"), py::arg("weights"), py::arg("seed"), py::arg("standardize") = false);

    py::class_<MarginalFactory>(m, "MarginalFactory")
        .def_static(
            "create",
            [](MarginalType type, py::object cfg_obj, py::object seed_obj) {
                unsigned int seed = seed_obj.is_none() ? fresh_seed() : seed_obj.cast<unsigned int>();
                MarginalConfig cfg = to_marginal_config(cfg_obj);
                auto up = MarginalFactory::create(type, std::move(cfg), seed);
                return std::shared_ptr<IMarginalDistribution>(up.release());
            },
            py::arg("type"), py::arg("config"), py::arg("seed") = py::none(),
            "Create a marginal distribution from a type and a typed config."
        )
        .def_static("create_gaussian", [](GaussianMarginalCfg cfg, py::object seed_obj) {
            unsigned int seed = seed_obj.is_none() ? fresh_seed() : seed_obj.cast<unsigned int>();
            auto up = MarginalFactory::create(MarginalType::GAUSSIAN, MarginalConfig{cfg}, seed);
            return std::shared_ptr<IMarginalDistribution>(up.release());
        }, py::arg("cfg"), py::arg("seed") = py::none())
        .def_static("create_split_gaussian", [](SplitGaussianMarginalCfg cfg, py::object seed_obj) {
            unsigned int seed = seed_obj.is_none() ? fresh_seed() : seed_obj.cast<unsigned int>();
            auto up = MarginalFactory::create(MarginalType::HALF_GAUSSIAN, MarginalConfig{cfg}, seed);
            return std::shared_ptr<IMarginalDistribution>(up.release());
        }, py::arg("cfg"), py::arg("seed") = py::none())
        .def_static("create_flat", [](FlatMarginalCfg cfg, py::object seed_obj) {
            unsigned int seed = seed_obj.is_none() ? fresh_seed() : seed_obj.cast<unsigned int>();
            auto up = MarginalFactory::create(MarginalType::FLAT, MarginalConfig{cfg}, seed);
            return std::shared_ptr<IMarginalDistribution>(up.release());
        }, py::arg("cfg"), py::arg("seed") = py::none())
        .def_static("create_likelihood", [](LikelihoodMarginalCfg cfg, py::object seed_obj) {
            unsigned int seed = seed_obj.is_none() ? fresh_seed() : seed_obj.cast<unsigned int>();
            auto up = MarginalFactory::create(MarginalType::LIKELIHOOD, MarginalConfig{cfg}, seed);
            return std::shared_ptr<IMarginalDistribution>(up.release());
        }, py::arg("cfg"), py::arg("seed") = py::none());
}

void init_marginal_config_factory(py::module_& m) {
    py::class_<MarginalConfigFactory>(m, "MarginalConfigFactory")
        .def(py::init<>())
        .def(
            "create",
            py::overload_cast<ParamId, MarginalType>(&MarginalConfigFactory::create),
            py::arg("pid"), py::arg("marginal"),
            "Create a MarginalConfig from a ParamId and a marginal type."
        )
        .def(
            "create",
            py::overload_cast<ExperimentObs, MarginalType>(&MarginalConfigFactory::create),
            py::arg("obs"), py::arg("marginal"),
            "Create a MarginalConfig from an ExperimentObs and a marginal type."
        );
}

void init_copulas(py::module_& m) {
    py::enum_<CopulaType>(m, "CopulaType")
        .value("GAUSSIAN", CopulaType::GAUSSIAN)
        .value("STUDENT_T", CopulaType::STUDENT_T)
        .export_values();

    py::class_<GaussianCopulaConfig, AbstractConfig>(m, "GaussianCopulaConfig")
        .def(py::init<>())
        .def(py::init([](decltype(GaussianCopulaConfig::R) R) {
            GaussianCopulaConfig c;
            c.R = std::move(R);
            return c;
        }), py::arg("R"))
        .def_readwrite("R", &GaussianCopulaConfig::R)
        .def("__repr__", [](const GaussianCopulaConfig&) {
            return std::string("GaussianCopulaConfig(R=...)");
        });

    py::class_<StudentTCopulaConfig, AbstractConfig>(m, "StudentTCopulaConfig")
        .def(py::init([]() {
            StudentTCopulaConfig c;
            c.nu = 4;
            return c;
        }))
        .def(py::init([](decltype(StudentTCopulaConfig::R) R, int nu) {
            StudentTCopulaConfig c;
            c.R = std::move(R);
            c.nu = nu;
            return c;
        }), py::arg("R"), py::arg("nu"))
        .def_readwrite("R", &StudentTCopulaConfig::R)
        .def_readwrite("nu", &StudentTCopulaConfig::nu)
        .def("__repr__", [](const StudentTCopulaConfig& c) {
            std::ostringstream oss;
            oss << "StudentTCopulaConfig(nu=" << c.nu << ", R=...)";
            return oss.str();
        });

    py::class_<ICopula, std::shared_ptr<ICopula>>(m, "ICopula")
        .def("sample_u",
             py::overload_cast<std::size_t>(&ICopula::sample_u),
             py::arg("n"))
        .def("sample_u",
             py::overload_cast<>(&ICopula::sample_u))
        .def("log_density", &ICopula::log_density, py::arg("u"));

    py::class_<GaussianCopula, ICopula, std::shared_ptr<GaussianCopula>>(m, "GaussianCopula")
        .def(py::init([](decltype(GaussianCopulaConfig::R) R, py::object seed_obj) {
            unsigned int seed = seed_obj.is_none() ? fresh_seed() : seed_obj.cast<unsigned int>();
            return std::make_shared<GaussianCopula>(seed, std::move(R));
        }), py::arg("R"), py::arg("seed") = py::none());

    py::class_<StudentTCopula, ICopula, std::shared_ptr<StudentTCopula>>(m, "StudentTCopula")
        .def(py::init([](decltype(StudentTCopulaConfig::R) R, int nu, py::object seed_obj) {
            unsigned int seed = seed_obj.is_none() ? fresh_seed() : seed_obj.cast<unsigned int>();
            return std::make_shared<StudentTCopula>(seed, std::move(R), nu);
        }), py::arg("R"), py::arg("nu"), py::arg("seed") = py::none());

    py::class_<CopulaFactory>(m, "CopulaFactory")
        .def_static(
            "create",
            [](CopulaType type, py::object config_obj, py::object seed_obj) {
                unsigned int seed = seed_obj.is_none() ? fresh_seed() : seed_obj.cast<unsigned int>();
                CopulaConfig cfg = to_copula_config(config_obj);

                if (py::isinstance<GaussianCopulaConfig>(config_obj) && type != CopulaType::GAUSSIAN) {
                    throw py::value_error("CopulaType/config mismatch: GaussianCopulaConfig requires CopulaType.GAUSSIAN.");
                }
                if (py::isinstance<StudentTCopulaConfig>(config_obj) && type != CopulaType::STUDENT_T) {
                    throw py::value_error("CopulaType/config mismatch: StudentTCopulaConfig requires CopulaType.STUDENT_T.");
                }

                auto up = CopulaFactory::create(type, std::move(cfg), seed);
                return std::shared_ptr<ICopula>(up.release());
            },
            py::arg("type"), py::arg("config"), py::arg("seed") = py::none(),
            "Create a copula from a type and a typed config."
        )
        .def_static("create_gaussian", [](GaussianCopulaConfig cfg, py::object seed_obj) {
            unsigned int seed = seed_obj.is_none() ? fresh_seed() : seed_obj.cast<unsigned int>();
            auto up = CopulaFactory::create(CopulaType::GAUSSIAN, CopulaConfig{cfg}, seed);
            return std::shared_ptr<ICopula>(up.release());
        }, py::arg("cfg"), py::arg("seed") = py::none())
        .def_static("create_student_t", [](StudentTCopulaConfig cfg, py::object seed_obj) {
            unsigned int seed = seed_obj.is_none() ? fresh_seed() : seed_obj.cast<unsigned int>();
            auto up = CopulaFactory::create(CopulaType::STUDENT_T, CopulaConfig{cfg}, seed);
            return std::shared_ptr<ICopula>(up.release());
        }, py::arg("cfg"), py::arg("seed") = py::none());
}

void init_joint_distribution(py::module_& m) {
    py::class_<JointDistribution, std::shared_ptr<JointDistribution>>(m, "JointDistribution")
        .def("sample",
             py::overload_cast<std::size_t>(&JointDistribution::sample, py::const_),
             py::arg("n"),
             "Sample n points from the joint distribution.")
        .def("sample",
             py::overload_cast<>(&JointDistribution::sample, py::const_),
             "Sample one point from the joint distribution.")
        .def("logpdf", &JointDistribution::logpdf, py::arg("x"))
        .def("curvature", &JointDistribution::curvature, py::arg("x"))
        .def("dim", &JointDistribution::dim)
        .def("get_stds", &JointDistribution::get_stds)
        .def_static(
            "create",
            [](const std::vector<MarginalType>& m_types,
               const std::vector<py::object>& m_cfgs,
               CopulaType c_type,
               py::object c_cfg,
               py::object seed_obj) {
                if (m_types.size() != m_cfgs.size()) {
                    throw py::value_error("marginal_types and marginal_configs must have the same length.");
                }
                if (m_types.empty()) {
                    throw py::value_error("At least one marginal is required.");
                }

                unsigned int base_seed = seed_obj.is_none() ? fresh_seed() : seed_obj.cast<unsigned int>();

                std::vector<std::unique_ptr<IMarginalDistribution>> marginals;
                marginals.reserve(m_types.size());
                for (std::size_t i = 0; i < m_types.size(); ++i) {
                    MarginalConfig cfg = to_marginal_config(m_cfgs[i]);
                    auto up = MarginalFactory::create(
                        m_types[i],
                        std::move(cfg),
                        base_seed + static_cast<unsigned int>(i)
                    );
                    marginals.push_back(std::move(up));
                }

                CopulaConfig ccfg = to_copula_config(c_cfg);
                auto copula = CopulaFactory::create(c_type, std::move(ccfg), base_seed + 1000003u);
                return std::make_shared<JointDistribution>(std::move(marginals), std::move(copula));
            },
            py::arg("marginal_types"),
            py::arg("marginal_configs"),
            py::arg("copula_type"),
            py::arg("copula_config"),
            py::arg("seed") = py::none(),
            "Create a JointDistribution from marginal configs and a copula config."
        )
        .def_static(
            "create_with_seeds",
            [](const std::vector<MarginalType>& m_types,
               const std::vector<py::object>& m_cfgs,
               const std::vector<unsigned int>& m_seeds,
               CopulaType c_type,
               py::object c_cfg,
               unsigned int copula_seed) {
                if (m_types.size() != m_cfgs.size() || m_types.size() != m_seeds.size()) {
                    throw py::value_error("marginal_types, marginal_configs and marginal_seeds must have the same length.");
                }
                if (m_types.empty()) {
                    throw py::value_error("At least one marginal is required.");
                }

                std::vector<std::unique_ptr<IMarginalDistribution>> marginals;
                marginals.reserve(m_types.size());
                for (std::size_t i = 0; i < m_types.size(); ++i) {
                    MarginalConfig cfg = to_marginal_config(m_cfgs[i]);
                    auto up = MarginalFactory::create(m_types[i], std::move(cfg), m_seeds[i]);
                    marginals.push_back(std::move(up));
                }

                CopulaConfig ccfg = to_copula_config(c_cfg);
                auto copula = CopulaFactory::create(c_type, std::move(ccfg), copula_seed);
                return std::make_shared<JointDistribution>(std::move(marginals), std::move(copula));
            },
            py::arg("marginal_types"),
            py::arg("marginal_configs"),
            py::arg("marginal_seeds"),
            py::arg("copula_type"),
            py::arg("copula_config"),
            py::arg("copula_seed"),
            "Create a JointDistribution with explicit seeds."
        );
}

void init_likelihood_and_profiling(py::module_& m) {
    py::class_<ILikelihood, std::shared_ptr<ILikelihood>>(m, "ILikelihood")
        .def("nll", &ILikelihood::nll, py::arg("theta"))
        .def("dim", &ILikelihood::dim);

    py::class_<BaseLikelihood, ILikelihood, std::shared_ptr<BaseLikelihood>>(m, "BaseLikelihood")
        .def("enable_debug_trace", &BaseLikelihood::enable_debug_trace, py::arg("max_evals") = 25)
        .def("disable_debug_trace", &BaseLikelihood::disable_debug_trace);

    py::class_<FitResult>(m, "FitResult")
        .def(py::init<>())
        .def_readwrite("p_hat", &FitResult::p_hat)
        .def_readwrite("eta_hat", &FitResult::eta_hat)
        .def_readwrite("p_hat_std", &FitResult::p_hat_std)
        .def_readwrite("p_hat_correlations", &FitResult::p_hat_correlations)
        .def_readwrite("ell_hat", &FitResult::ell_hat);

    py::class_<ProfileRequest>(m, "ProfileRequest")
        .def(py::init<>())
        .def_readwrite("free_params", &ProfileRequest::free_params)
        .def_readwrite("fixed_params", &ProfileRequest::fixed_params)
        .def_readwrite("start", &ProfileRequest::start);

    py::class_<ProfileResult>(m, "ProfileResult")
        .def(py::init<>())
        .def_readwrite("nll_hat", &ProfileResult::nll_hat)
        .def_readwrite("theta_hat", &ProfileResult::theta_hat)
        .def_readwrite("converged", &ProfileResult::converged);

#if STAT_BIND_HAS_PROFILER
    if (!py::hasattr(m, "IFitBackend")) {
        py::class_<fit_app::IFitBackend, std::shared_ptr<fit_app::IFitBackend>>(m, "IFitBackend");
    }

    py::class_<Profiler>(m, "Profiler")
        .def(py::init<std::shared_ptr<fit_app::IFitBackend>>(), py::arg("minimizer"))
        .def("profile", &Profiler::profile, py::arg("base"), py::arg("request"));
#endif

#if STAT_BIND_HAS_GAUSSIAN_CONSTRAINTS
    py::class_<WithGaussianConstraints, ILikelihood, std::shared_ptr<WithGaussianConstraints>>(m, "WithGaussianConstraints")
        .def(py::init<std::shared_ptr<ILikelihood>, std::shared_ptr<JointDistribution>, std::vector<std::size_t>>(),
             py::arg("base"), py::arg("constraints_dist"), py::arg("constrained_params"));
#endif
}

void init_fit_and_contours(py::module_& m) {
    py::enum_<ProfilingMethod>(m, "ProfilingMethod")
        .value("SLICE", ProfilingMethod::SLICE)
        .value("FREE_PROJECTION", ProfilingMethod::FREE_PROJECTION)
        .value("PRIOR_CONSTRAINED_PROJECTION", ProfilingMethod::PRIOR_CONSTRAINED_PROJECTION)
        .export_values();

    py::enum_<ContourAlgorithm>(m, "ContourAlgorithm")
        .value("AMS", ContourAlgorithm::AMS)
        .value("MINUIT", ContourAlgorithm::MINUIT)
        .export_values();

    py::enum_<ProfileBackend>(m, "ProfileBackend")
        .value("MINUIT", ProfileBackend::MINUIT)
        .value("LAPLACE_NUISANCE", ProfileBackend::LAPLACE_NUISANCE)
        .export_values();

    py::class_<MLFitOptions>(m, "MLFitOptions")
        .def(py::init<>())
        .def_readwrite("run_hesse", &MLFitOptions::run_hesse)
        .def_readwrite("request_minos", &MLFitOptions::request_minos)
        .def_readwrite("verbose", &MLFitOptions::verbose)
        .def_readwrite("strategy", &MLFitOptions::strategy)
        .def_readwrite("max_fcn", &MLFitOptions::max_fcn)
        .def_readwrite("tolerance", &MLFitOptions::tolerance)
        .def_readwrite("allow_profile_hessian_fallback", &MLFitOptions::allow_profile_hessian_fallback)
        .def_readwrite("profile_hessian_step_scale", &MLFitOptions::profile_hessian_step_scale)
        .def_readwrite("profile_hessian_eig_floor_rel", &MLFitOptions::profile_hessian_eig_floor_rel)
        .def_readwrite("trace_first_evals", &MLFitOptions::trace_first_evals)
        .def_readwrite("trace_max_evals", &MLFitOptions::trace_max_evals);

    py::class_<ContourOptions>(m, "ContourOptions")
        .def(py::init<>())
        .def_readwrite("profiling_method", &ContourOptions::profiling_method)
        .def_readwrite("profile_backend", &ContourOptions::profile_backend)
        .def_readwrite("primary_contour_method", &ContourOptions::primary_contour_method)
        .def_readwrite("fallback_contour_method", &ContourOptions::fallback_contour_method)
        .def_readwrite("resolution", &ContourOptions::resolution);
        // on_progress n'est volontairement pas exposé ici : sa signature dépend
        // de ContourObserver.h. Ajoute un binding dédié si tu veux un callback Python.

    py::class_<Contour>(m, "Contour")
        .def_property_readonly("paths", [](const Contour& contour) {
            return std::vector<Path>(contour.paths.begin(), contour.paths.end());
        })
        .def_readonly("level", &Contour::level)
        .def_readonly("success", &Contour::success)
        .def("__repr__", [](const Contour& contour) {
            std::ostringstream os;
            os << "<Contour success=" << (contour.success ? "true" : "false")
               << " level=" << contour.level
               << " paths=" << contour.paths.size() << ">";
            return os.str();
        });
}

void init_statistic_data_structures(py::module_& m) {
    py::enum_<StatisticLikelihoodMode>(m, "StatisticLikelihoodMode")
        .value("PROFILED_NUISANCE", StatisticLikelihoodMode::PROFILED_NUISANCE)
        .value("CHI2_MC_COVARIANCE", StatisticLikelihoodMode::CHI2_MC_COVARIANCE)
        .export_values();

    py::class_<MCConfig>(m, "MCConfig")
        .def(py::init<>())
        .def_readwrite("draws", &MCConfig::draws)
        .def_readwrite("skew_abs_threshold", &MCConfig::skew_abs_threshold)
        .def_readwrite("covariance_ridge_rel", &MCConfig::covariance_ridge_rel)
        .def_readwrite("covariance_ridge_abs", &MCConfig::covariance_ridge_abs)
        .def_readwrite("retry_failed_predictions", &MCConfig::retry_failed_predictions)
        .def_readwrite("max_prediction_failures", &MCConfig::max_prediction_failures);

    py::class_<MCObservableCovariance>(m, "MCObservableCovariance")
        .def(py::init<>())
        .def_readwrite("ids", &MCObservableCovariance::ids)
        .def_readwrite("mean", &MCObservableCovariance::mean)
        .def_readwrite("covariance", &MCObservableCovariance::covariance)
        .def_readwrite("covariance_inv", &MCObservableCovariance::covariance_inv);

    m.def("covariance_from_obs_samples", &covariance_from_obs_samples,
          py::arg("samples"), py::arg("ids"),
          py::arg("ridge_rel") = 1e-8, py::arg("ridge_abs") = 1e-12);
    m.def("covariance_ids_from_first_sample", &covariance_ids_from_first_sample,
          py::arg("samples"));

    py::class_<GaussianSummary>(m, "GaussianSummary")
        .def(py::init<>())
        .def_readwrite("id", &GaussianSummary::id)
        .def_readwrite("mu", &GaussianSummary::mu)
        .def_readwrite("sigma", &GaussianSummary::sigma)
        .def_readwrite("sigma_p", &GaussianSummary::sigma_p)
        .def_readwrite("sigma_m", &GaussianSummary::sigma_m)
        .def_readwrite("mode", &GaussianSummary::mode)
        .def_readwrite("skew", &GaussianSummary::skew)
        .def_readwrite("symmetric", &GaussianSummary::symmetric)
        .def("__repr__", [](const GaussianSummary& gs) {
            std::ostringstream oss;
            oss << gs;
            return "<GaussianSummary " + oss.str() + ">";
        })
        .def("__str__", [](const GaussianSummary& gs) {
            std::ostringstream oss;
            oss << gs;
            return oss.str();
        });

    py::class_<MCRealization>(m, "MCRealization")
        .def(py::init<>())
        .def_readwrite("sampled_obss", &MCRealization::sampled_obss)
        .def_readwrite("sampled_params", &MCRealization::sampled_params);

    py::class_<MCResult>(m, "MCResult")
        .def(py::init<>())
        .def_readwrite("mc_real", &MCResult::mc_real)
        .def_readwrite("summary", &MCResult::summary)
        .def_readwrite("covariance", &MCResult::covariance);

    py::class_<AdvancedStatisticConfig>(m, "AdvancedStatisticConfig", R"pbdoc(
Advanced statistical configuration.

This object groups fit, nuisance-pruning, covariance and likelihood-backend
options that are useful for expert workflows but too specialized for the basic
``StatisticConfig`` surface.
)pbdoc")
        .def(py::init<>())
        .def_readwrite("override_nuisance_marginals", &AdvancedStatisticConfig::override_nuisance_marginals)
        .def_readwrite("override_exp_data_marginals", &AdvancedStatisticConfig::override_exp_data_marginals)
        .def_readwrite("nuisance_copula_type", &AdvancedStatisticConfig::nuisance_copula_type)
        .def_readwrite("exp_data_copula_type", &AdvancedStatisticConfig::exp_data_copula_type)
        .def_readwrite("MLE_max_iter", &AdvancedStatisticConfig::MLE_max_iter)
        .def_readwrite("MLE_tol", &AdvancedStatisticConfig::MLE_tol)
        .def_readwrite("MLE_strategy", &AdvancedStatisticConfig::MLE_strategy)
        .def_readwrite("MLE_run_hesse", &AdvancedStatisticConfig::MLE_run_hesse)
        .def_readwrite("MLE_request_minos", &AdvancedStatisticConfig::MLE_request_minos)
        .def_readwrite("MLE_verbose", &AdvancedStatisticConfig::MLE_verbose)
        .def_readwrite("nuisance_relevance_cutoff", &AdvancedStatisticConfig::nuisance_relevance_cutoff)
        .def_readwrite("nuisance_sensitivity_pruning", &AdvancedStatisticConfig::nuisance_sensitivity_pruning)
        .def_readwrite("nuisance_sensitivity_probe_sigmas", &AdvancedStatisticConfig::nuisance_sensitivity_probe_sigmas)
        .def_readwrite("nuisance_sensitivity_rel_cutoff", &AdvancedStatisticConfig::nuisance_sensitivity_rel_cutoff)
        .def_readwrite("nuisance_sensitivity_abs_cutoff", &AdvancedStatisticConfig::nuisance_sensitivity_abs_cutoff)
        .def_readwrite("nuisance_sensitivity_scale_floor", &AdvancedStatisticConfig::nuisance_sensitivity_scale_floor)
        .def_readwrite("MLE_trace_first_evals", &AdvancedStatisticConfig::MLE_trace_first_evals)
        .def_readwrite("MLE_trace_max_evals", &AdvancedStatisticConfig::MLE_trace_max_evals)
        .def_readwrite("MLE_allow_profile_hessian_fallback", &AdvancedStatisticConfig::MLE_allow_profile_hessian_fallback)
        .def_readwrite("MLE_profile_hessian_step_scale", &AdvancedStatisticConfig::MLE_profile_hessian_step_scale)
        .def_readwrite("MLE_profile_hessian_eig_floor_rel", &AdvancedStatisticConfig::MLE_profile_hessian_eig_floor_rel)
        .def_readwrite("likelihood_mode", &AdvancedStatisticConfig::likelihood_mode)
        .def_readwrite("chi2_covariance_ridge_rel", &AdvancedStatisticConfig::chi2_covariance_ridge_rel)
        .def_readwrite("chi2_covariance_ridge_abs", &AdvancedStatisticConfig::chi2_covariance_ridge_abs)
        .def_readwrite("MC_force_decay_threads_to_one", &AdvancedStatisticConfig::MC_force_decay_threads_to_one)
        .def_readwrite("MC_forced_decay_threads", &AdvancedStatisticConfig::MC_forced_decay_threads)
        .def_readwrite("nuisance_sensitivity_contexts", &AdvancedStatisticConfig::nuisance_sensitivity_contexts)
        .def_readwrite("nuisance_sensitivity_context_sigma", &AdvancedStatisticConfig::nuisance_sensitivity_context_sigma)
        .def_readwrite("nuisance_sensitivity_seed", &AdvancedStatisticConfig::nuisance_sensitivity_seed)
        .def_readwrite("nuisance_sensitivity_keep_on_failure", &AdvancedStatisticConfig::nuisance_sensitivity_keep_on_failure);


    py::class_<StatisticProgressEvent>(m, "StatisticProgressEvent")
        .def(py::init<>())
        .def_readwrite("phase", &StatisticProgressEvent::phase)
        .def_readwrite("message", &StatisticProgressEvent::message)
        .def_readwrite("completed", &StatisticProgressEvent::completed)
        .def_readwrite("total", &StatisticProgressEvent::total)
        .def_readwrite("attempts", &StatisticProgressEvent::attempts)
        .def_readwrite("failures", &StatisticProgressEvent::failures)
        .def_readwrite("fraction", &StatisticProgressEvent::fraction)
        .def_readwrite("elapsed_seconds", &StatisticProgressEvent::elapsed_seconds)
        .def_readwrite("eta_seconds", &StatisticProgressEvent::eta_seconds)
        .def_readwrite("finished", &StatisticProgressEvent::finished)
        .def_readwrite("sequence", &StatisticProgressEvent::sequence);

    py::class_<StatisticProgressMonitor, std::shared_ptr<StatisticProgressMonitor>>(m, "StatisticProgressMonitor")
        .def(py::init<>())
        .def("reset", &StatisticProgressMonitor::reset,
             py::arg("phase") = "preparing",
             py::arg("message") = "Preparing statistic workflow")
        .def("snapshot", &StatisticProgressMonitor::snapshot)
        .def("set_progress", [](StatisticProgressMonitor& monitor,
                                const std::string& phase,
                                const std::string& message,
                                double fraction,
                                std::size_t completed,
                                std::size_t total,
                                double eta_seconds,
                                bool finished) {
            StatisticProgressEvent event;
            event.phase = phase;
            event.message = message;
            event.fraction = fraction;
            event.completed = completed;
            event.total = total;
            event.eta_seconds = eta_seconds;
            event.finished = finished;
            monitor.update(std::move(event));
        },
        py::arg("phase"), py::arg("message"), py::arg("fraction"),
        py::arg("completed") = 0, py::arg("total") = 0,
        py::arg("eta_seconds") = -1.0, py::arg("finished") = false);

    py::class_<StatisticConfig>(m, "StatisticConfig", R"pbdoc(
Basic statistical configuration.

Keep common runtime controls here: MC draw count, skew threshold, print toggles
and output options.  Expert knobs live under ``advanced``.
)pbdoc")
        .def(py::init<>())
        .def_readwrite("MC_draws", &StatisticConfig::MC_draws)
        .def_readwrite("MC_threads", &StatisticConfig::MC_threads)
        .def_readwrite("MC_seed", &StatisticConfig::MC_seed)
        .def_readwrite("skew_abs_threshold", &StatisticConfig::skew_abs_threshold)
        .def_readwrite("print_mc_progress", &StatisticConfig::print_mc_progress)
        .def_readwrite("print_chi2_pipeline_progress", &StatisticConfig::print_chi2_pipeline_progress)
        .def_readwrite("print_mc_config", &StatisticConfig::print_mc_config)
        .def_readwrite("print_fit_summary", &StatisticConfig::print_fit_summary)
        .def_readwrite("print_scan_summary", &StatisticConfig::print_scan_summary)
        .def_readwrite("print_cache_summary", &StatisticConfig::print_cache_summary)
        .def_readwrite("print_debug", &StatisticConfig::print_debug)
        .def_readwrite("write_mc_samples_csv", &StatisticConfig::write_mc_samples_csv)
        .def_readwrite("mc_samples_csv_path", &StatisticConfig::mc_samples_csv_path)
        .def_readwrite("mc_progress_probe_draws", &StatisticConfig::mc_progress_probe_draws)
        .def_readwrite("mc_progress_update_every", &StatisticConfig::mc_progress_update_every)
        .def_readwrite("progress_monitor", &StatisticConfig::progress_monitor)
        .def_readwrite("fit_parameter_bounds", &StatisticConfig::fit_parameter_bounds)
        .def_readwrite("fit_parameter_offsets", &StatisticConfig::fit_parameter_offsets)
        .def_readwrite("advanced", &StatisticConfig::advanced);


    py::class_<FitResultWithMaps>(m, "FitResultWithMaps")
        .def(py::init<>())
        .def_readwrite("fit_ok", &FitResultWithMaps::fit_ok)
        .def_readwrite("p_hat", &FitResultWithMaps::p_hat)
        .def_readwrite("eta_hat", &FitResultWithMaps::eta_hat)
        .def_readwrite("p_hat_std", &FitResultWithMaps::p_hat_std)
        .def_readwrite("p_correlations", &FitResultWithMaps::p_correlations)
        .def_readwrite("ell_hat", &FitResultWithMaps::ell_hat);

    py::class_<LikelihoodScanPoint>(m, "LikelihoodScanPoint")
        .def(py::init<>())
        .def_readwrite("x", &LikelihoodScanPoint::x)
        .def_readwrite("y", &LikelihoodScanPoint::y)
        .def_readwrite("nll", &LikelihoodScanPoint::nll)
        .def_readwrite("delta_nll", &LikelihoodScanPoint::delta_nll);

    py::class_<LikelihoodScanGrid>(m, "LikelihoodScanGrid")
        .def(py::init<>())
        .def_readwrite("x_param", &LikelihoodScanGrid::x_param)
        .def_readwrite("y_param", &LikelihoodScanGrid::y_param)
        .def_readwrite("x_center", &LikelihoodScanGrid::x_center)
        .def_readwrite("y_center", &LikelihoodScanGrid::y_center)
        .def_readwrite("nx", &LikelihoodScanGrid::nx)
        .def_readwrite("ny", &LikelihoodScanGrid::ny)
        .def_readwrite("points", &LikelihoodScanGrid::points);
}

void init_statistic_interface(py::module_& m) {
    py::class_<StatisticInterface, std::shared_ptr<StatisticInterface>>(m, "StatisticInterface")
        .def(py::init<StatisticConfig, std::shared_ptr<ObservableInterface>>(),
             py::arg("config"), py::arg("observable_interface"))
        .def("select_experiment", &StatisticInterface::select_experiment, py::arg("experiment"))
        .def("select_experiments", &StatisticInterface::select_experiments, py::arg("experiments"))
        .def("select_experiments_all", &StatisticInterface::select_experiments_all)
        .def("has_experiment_selection", &StatisticInterface::has_experiment_selection)
        .def("selected_experiments", &StatisticInterface::selected_experiments)
        .def("compute_uncertainties", &StatisticInterface::compute_uncertainties,
             py::call_guard<py::gil_scoped_release>())
        .def("compute_uncertainties_and_sampling", &StatisticInterface::compute_uncertainties_and_sampling,
             py::call_guard<py::gil_scoped_release>())
        .def("compute_MLE", &StatisticInterface::compute_MLE, py::arg("p_specs"),
             py::call_guard<py::gil_scoped_release>())
        .def("compute_confidence_contour",
             &StatisticInterface::compute_confidence_contour,
             py::arg("p1"),
             py::arg("p2"),
             py::arg("z"),
             py::arg("bounds"),
             py::arg("options") = ContourOptions{},
             py::call_guard<py::gil_scoped_release>())
        .def("reload_nuisance_specs", &StatisticInterface::reload_nuisance_specs)
        .def("set_nuisance_user_file", &StatisticInterface::set_nuisance_user_file,
             py::arg("user_yaml_path"))
        .def("clear_nuisance_user_file", &StatisticInterface::clear_nuisance_user_file)
        .def("prepare_likelihood_for_scan", &StatisticInterface::prepare_likelihood_for_scan,
             py::arg("p_specs"))
        .def("set_manual_scan_point", &StatisticInterface::set_manual_scan_point,
             py::arg("p_hat"), py::arg("eta_hat"))
        .def("scan_likelihood_around_current_point",
             &StatisticInterface::scan_likelihood_around_current_point,
             py::arg("p1"),
             py::arg("p2"),
             py::arg("x_half_width"),
             py::arg("y_half_width"),
             py::arg("nx"),
             py::arg("ny"))
        .def("save_likelihood_scan_csv", &StatisticInterface::save_likelihood_scan_csv,
             py::arg("path"), py::arg("grid"))
        .def("update_cache", &StatisticInterface::update_cache,
             py::arg("p_specs") = std::vector<ParamId>{})
        .def("get_all_obss_deps", &StatisticInterface::get_all_obss_deps)
        .def("get_active_observable_dependencies", &StatisticInterface::get_active_observable_dependencies,
             R"pbdoc(Return the parameter dependencies currently visible to the Statistic manager.)pbdoc")
        .def("get_p_specs", &StatisticInterface::get_p_specs,
             py::arg("p_specs"))
        .def("get_all_correlations", &StatisticInterface::get_all_correlations)
        .def("get_all_obs_correlations", &StatisticInterface::get_all_obs_correlations)
        .def("get_obs_exp", &StatisticInterface::get_obs_exp)
        .def("print_cache", &StatisticInterface::print_cache);
}

} // namespace

void init_statistic(py::module_& m) {
    // init_abstract_config(m);
    init_marginals(m);
    init_marginal_config_factory(m);
    init_copulas(m);
    init_joint_distribution(m);
    init_likelihood_and_profiling(m);
    init_fit_and_contours(m);
    init_statistic_data_structures(m);
    init_statistic_interface(m);
}

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
        .export_values();

    py::enum_<ContourAlgorithm>(m, "ContourAlgorithm")
        .value("MINUIT", ContourAlgorithm::MINUIT)
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
        .def_readwrite("primary_contour_method", &ContourOptions::primary_contour_method)
        .def_readwrite("fallback_contour_method", &ContourOptions::fallback_contour_method)
        .def_readwrite("resolution", &ContourOptions::resolution);
        // on_progress n'est volontairement pas exposé ici : sa signature dépend
        // de ContourObserver.h. Ajoute un binding dédié si tu veux un callback Python.

    // Enregistrement minimal pour que StatisticInterface.compute_confidence_contour()
    // puisse retourner un objet Contour. Ajoute ici les .def_readwrite(...) si ton
    // struct Contour expose des champs publics.
    py::class_<Contour>(m, "Contour")
        .def("__repr__", [](const Contour&) { return std::string("<Contour>"); });
}

void init_statistic_data_structures(py::module_& m) {
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
        .def_readwrite("summary", &MCResult::summary);

    py::class_<StatisticConfig>(m, "StatisticConfig")
        .def(py::init<>())
        .def_readwrite("override_nuisance_marginals", &StatisticConfig::override_nuisance_marginals)
        .def_readwrite("override_exp_data_marginals", &StatisticConfig::override_exp_data_marginals)
        .def_readwrite("nuisance_copula_type", &StatisticConfig::nuisance_copula_type)
        .def_readwrite("exp_data_copula_type", &StatisticConfig::exp_data_copula_type)
        .def_readwrite("MC_draws", &StatisticConfig::MC_draws)
        .def_readwrite("skew_abs_threshold", &StatisticConfig::skew_abs_threshold)
        .def_readwrite("MLE_max_iter", &StatisticConfig::MLE_max_iter)
        .def_readwrite("MLE_tol", &StatisticConfig::MLE_tol)
        .def_readwrite("MLE_strategy", &StatisticConfig::MLE_strategy)
        .def_readwrite("MLE_run_hesse", &StatisticConfig::MLE_run_hesse)
        .def_readwrite("MLE_request_minos", &StatisticConfig::MLE_request_minos)
        .def_readwrite("MLE_verbose", &StatisticConfig::MLE_verbose)
        .def_readwrite("nuisance_relevance_cutoff", &StatisticConfig::nuisance_relevance_cutoff)
        .def_readwrite("nuisance_sensitivity_pruning", &StatisticConfig::nuisance_sensitivity_pruning)
        .def_readwrite("nuisance_sensitivity_probe_sigmas", &StatisticConfig::nuisance_sensitivity_probe_sigmas)
        .def_readwrite("nuisance_sensitivity_rel_cutoff", &StatisticConfig::nuisance_sensitivity_rel_cutoff)
        .def_readwrite("nuisance_sensitivity_abs_cutoff", &StatisticConfig::nuisance_sensitivity_abs_cutoff)
        .def_readwrite("nuisance_sensitivity_scale_floor", &StatisticConfig::nuisance_sensitivity_scale_floor)
        .def_readwrite("MLE_trace_first_evals", &StatisticConfig::MLE_trace_first_evals)
        .def_readwrite("MLE_trace_max_evals", &StatisticConfig::MLE_trace_max_evals)
        .def_readwrite("MLE_allow_profile_hessian_fallback", &StatisticConfig::MLE_allow_profile_hessian_fallback)
        .def_readwrite("MLE_profile_hessian_step_scale", &StatisticConfig::MLE_profile_hessian_step_scale)
        .def_readwrite("MLE_profile_hessian_eig_floor_rel", &StatisticConfig::MLE_profile_hessian_eig_floor_rel);

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
        .def("compute_uncertainties", &StatisticInterface::compute_uncertainties)
        .def("compute_uncertainties_and_sampling", &StatisticInterface::compute_uncertainties_and_sampling)
        .def("compute_MLE", &StatisticInterface::compute_MLE, py::arg("p_specs"))
        .def("compute_confidence_contour",
             &StatisticInterface::compute_confidence_contour,
             py::arg("p1"),
             py::arg("p2"),
             py::arg("z"),
             py::arg("bounds"),
             py::arg("options") = ContourOptions{});
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

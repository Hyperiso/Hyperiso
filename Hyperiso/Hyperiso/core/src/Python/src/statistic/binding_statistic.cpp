#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "StatisticInterface.h"

#include <sstream>
#include <random>
#include <memory>
#include "AbstractConfig.h"

#include "MarginalType.h"
#include "IMarginalDistribution.h"

#include "GaussianMarginal.h"
#include "SplitGaussianMarginal.h"
#include "FlatMarginal.h"
#include "LikelihoodMarginal.h"
#include "CopulaType.h"
#include "ICopula.h"
#include "GaussianCopula.h"
#include "StudentTCopula.h"
#include "CopulaFactory.h"
#include "MarginalFactory.h"
#include "JointDistribution.h"

namespace py = pybind11;

static unsigned int fresh_seed() {
    return std::random_device{}();
}

void init_marginals(py::module_ &m) {

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
        .def("__repr__", [](const FlatMarginalCfg &c) {
            std::ostringstream oss;
            oss << "FlatMarginalCfg(a=" << c.a << ", b=" << c.b << ")";
            return oss.str();
        });

    py::class_<GaussianMarginalCfg, AbstractConfig>(m, "GaussianMarginalCfg")
        .def(py::init<double, double>(), py::arg("mu"), py::arg("sigma"))
        .def_readwrite("mu", &GaussianMarginalCfg::mu)
        .def_readwrite("sigma", &GaussianMarginalCfg::sigma)
        .def("__repr__", [](const GaussianMarginalCfg &c) {
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
        .def("__repr__", [](const SplitGaussianMarginalCfg &c) {
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
        .def("__repr__", [](const LikelihoodMarginalCfg &c) {
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


    py::class_<DistributionFactory>(m, "DistributionFactory")
        .def_static(
            "create",
            [](MarginalType type, MarginalConfig cfg, py::object seed_obj) {
                unsigned int seed = seed_obj.is_none() ? fresh_seed() : seed_obj.cast<unsigned int>();
                auto up = DistributionFactory::create(type, std::move(cfg), seed);
                return std::shared_ptr<IMarginalDistribution>(up.release());
            },
            py::arg("type"),
            py::arg("config"),
            py::arg("seed") = py::none(),
            "Create a marginal distribution from a type and a config.\n"
            "config must be one of: FlatMarginalCfg, GaussianMarginalCfg, SplitGaussianMarginalCfg, LikelihoodMarginalCfg."
        )

        .def_static("create_gaussian", [](GaussianMarginalCfg cfg, py::object seed_obj) {
            unsigned int seed = seed_obj.is_none() ? fresh_seed() : seed_obj.cast<unsigned int>();
            auto up = DistributionFactory::create(MarginalType::GAUSSIAN, MarginalConfig{cfg}, seed);
            return std::shared_ptr<IMarginalDistribution>(up.release());
        }, py::arg("cfg"), py::arg("seed") = py::none())

        .def_static("create_split_gaussian", [](SplitGaussianMarginalCfg cfg, py::object seed_obj) {
            unsigned int seed = seed_obj.is_none() ? fresh_seed() : seed_obj.cast<unsigned int>();
            auto up = DistributionFactory::create(MarginalType::HALF_GAUSSIAN, MarginalConfig{cfg}, seed);
            return std::shared_ptr<IMarginalDistribution>(up.release());
        }, py::arg("cfg"), py::arg("seed") = py::none())

        .def_static("create_flat", [](FlatMarginalCfg cfg, py::object seed_obj) {
            unsigned int seed = seed_obj.is_none() ? fresh_seed() : seed_obj.cast<unsigned int>();
            auto up = DistributionFactory::create(MarginalType::FLAT, MarginalConfig{cfg}, seed);
            return std::shared_ptr<IMarginalDistribution>(up.release());
        }, py::arg("cfg"), py::arg("seed") = py::none())

        .def_static("create_likelihood", [](LikelihoodMarginalCfg cfg, py::object seed_obj) {
            unsigned int seed = seed_obj.is_none() ? fresh_seed() : seed_obj.cast<unsigned int>();
            auto up = DistributionFactory::create(MarginalType::LIKELIHOOD, MarginalConfig{cfg}, seed);
            return std::shared_ptr<IMarginalDistribution>(up.release());
        }, py::arg("cfg"), py::arg("seed") = py::none());
}

namespace py = pybind11;

void init_marginal_config_factory(py::module_ &m) {
    py::class_<MarginalConfigFactory>(m, "MarginalConfigFactory")
        .def(py::init<>())

        .def(
            "create",
            py::overload_cast<ParamId, MarginalType>(&MarginalConfigFactory::create),
            py::arg("pid"),
            py::arg("marginal"),
            "Create a MarginalConfig from a ParamId + marginal type."
        )

        .def(
            "create",
            py::overload_cast<ObservableId, MarginalType>(&MarginalConfigFactory::create),
            py::arg("oid"),
            py::arg("marginal"),
            "Create a MarginalConfig from an ObservableId + marginal type."
        );
}

void init_copulas(py::module_ &m) {

    py::enum_<CopulaType>(m, "CopulaType")
        .value("GAUSSIAN", CopulaType::GAUSSIAN)
        .value("STUDENT_T", CopulaType::STUDENT_T)
        .export_values();

    py::class_<GaussianCopulaConfig, AbstractConfig>(m, "GaussianCopulaConfig")
        .def(py::init([]() {
            GaussianCopulaConfig c;
            return c;
        }))
        .def(py::init([](decltype(GaussianCopulaConfig::R) R) {
            GaussianCopulaConfig c;
            c.R = std::move(R);
            return c;
        }), py::arg("R"))
        .def_readwrite("R", &GaussianCopulaConfig::R)
        .def("__repr__", [](const GaussianCopulaConfig &) {
            return "GaussianCopulaConfig(R=...)";
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
        .def("__repr__", [](const StudentTCopulaConfig &c) {
            std::ostringstream oss;
            oss << "StudentTCopulaConfig(nu=" << c.nu << ", R=...)";
            return oss.str();
        });

    py::class_<ICopula, std::shared_ptr<ICopula>>(m, "ICopula")
        .def("sample_u",
             py::overload_cast<std::size_t>(&ICopula::sample_u),
             py::arg("n"),
             "Sample n points from the copula (returns List[List[float]] in Python).")
        .def("sample_u",
             py::overload_cast<>(&ICopula::sample_u),
             "Sample a single point (returns List[float] in Python).")
        .def("log_density", &ICopula::log_density, py::arg("u"));


    py::class_<GaussianCopula, ICopula, std::shared_ptr<GaussianCopula>>(m, "GaussianCopula")
        // Python-friendly: GaussianCopula(R, seed=None)
        .def(py::init([](decltype(GaussianCopulaConfig::R) R, py::object seed_obj) {
                unsigned int seed = seed_obj.is_none() ? fresh_seed() : seed_obj.cast<unsigned int>();
                return std::make_shared<GaussianCopula>(seed, std::move(R));
            }),
            py::arg("R"),
            py::arg("seed") = py::none()
        );

    py::class_<StudentTCopula, ICopula, std::shared_ptr<StudentTCopula>>(m, "StudentTCopula")
        // Python-friendly: StudentTCopula(R, nu, seed=None)
        .def(py::init([](decltype(StudentTCopulaConfig::R) R, int nu, py::object seed_obj) {
                unsigned int seed = seed_obj.is_none() ? fresh_seed() : seed_obj.cast<unsigned int>();
                return std::make_shared<StudentTCopula>(seed, std::move(R), nu);
            }),
            py::arg("R"),
            py::arg("nu"),
            py::arg("seed") = py::none()
        );

    py::class_<CopulaFactory>(m, "CopulaFactory")
        .def_static(
            "create",
            [](CopulaType type, py::object config_obj, py::object seed_obj) {
                unsigned int seed = seed_obj.is_none() ? fresh_seed() : seed_obj.cast<unsigned int>();

                CopulaConfig cfg;

                if (py::isinstance<GaussianCopulaConfig>(config_obj)) {
                    cfg = config_obj.cast<GaussianCopulaConfig>();
                    if (type != CopulaType::GAUSSIAN) {
                        throw py::value_error("CopulaType/config mismatch: got GaussianCopulaConfig but type != GAUSSIAN.");
                    }
                } else if (py::isinstance<StudentTCopulaConfig>(config_obj)) {
                    cfg = config_obj.cast<StudentTCopulaConfig>();
                    if (type != CopulaType::STUDENT_T) {
                        throw py::value_error("CopulaType/config mismatch: got StudentTCopulaConfig but type != STUDENT_T.");
                    }
                } else {
                    throw py::type_error("config must be GaussianCopulaConfig or StudentTCopulaConfig.");
                }

                auto up = CopulaFactory::create(type, std::move(cfg), seed);
                return std::shared_ptr<ICopula>(up.release());
            },
            py::arg("type"),
            py::arg("config"),
            py::arg("seed") = py::none(),
            "Create a copula from a type + config. Returns ICopula."
        )
        .def_static(
            "create_gaussian",
            [](GaussianCopulaConfig cfg, py::object seed_obj) {
                unsigned int seed = seed_obj.is_none() ? fresh_seed() : seed_obj.cast<unsigned int>();
                auto up = CopulaFactory::create(CopulaType::GAUSSIAN, CopulaConfig{cfg}, seed);
                return std::shared_ptr<ICopula>(up.release());
            },
            py::arg("cfg"),
            py::arg("seed") = py::none()
        )
        .def_static(
            "create_student_t",
            [](StudentTCopulaConfig cfg, py::object seed_obj) {
                unsigned int seed = seed_obj.is_none() ? fresh_seed() : seed_obj.cast<unsigned int>();
                auto up = CopulaFactory::create(CopulaType::STUDENT_T, CopulaConfig{cfg}, seed);
                return std::shared_ptr<ICopula>(up.release());
            },
            py::arg("cfg"),
            py::arg("seed") = py::none()
        );
}

static MarginalConfig to_marginal_config(const py::handle &obj) {
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
    throw py::type_error("Marginal config must be one of: FlatMarginalCfg, GaussianMarginalCfg, "
                         "SplitGaussianMarginalCfg, LikelihoodMarginalCfg.");
}

static CopulaConfig to_copula_config(const py::handle &obj) {
    if (py::isinstance<GaussianCopulaConfig>(obj)) {
        return obj.cast<GaussianCopulaConfig>();
    }
    if (py::isinstance<StudentTCopulaConfig>(obj)) {
        return obj.cast<StudentTCopulaConfig>();
    }
    throw py::type_error("Copula config must be one of: GaussianCopulaConfig, StudentTCopulaConfig.");
}

void init_joint_distribution(py::module_ &m) {
    py::class_<JointDistribution, std::shared_ptr<JointDistribution>>(m, "JointDistribution")
        .def("sample",
             py::overload_cast<std::size_t>(&JointDistribution::sample, py::const_),
             py::arg("n"),
             "Sample n points from the joint distribution. Returns List[List[float]].")
        .def("sample",
             py::overload_cast<>(&JointDistribution::sample, py::const_),
             "Sample a single point. Returns List[float].")

        .def("logpdf", &JointDistribution::logpdf, py::arg("x"))
        .def("dim", &JointDistribution::dim)

        .def_static(
            "create",
            [](const std::vector<MarginalType> &m_types,
               const std::vector<py::object> &m_cfgs,
               CopulaType c_type,
               py::object c_cfg,
               py::object seed_obj) {

                if (m_types.size() != m_cfgs.size()) {
                    throw py::value_error("marginal_types and marginal_configs must have the same length.");
                }
                if (m_types.empty()) {
                    throw py::value_error("At least one marginal is required.");
                }

                unsigned int base_seed = seed_obj.is_none() ? fresh_seed()
                                                            : seed_obj.cast<unsigned int>();

                std::vector<std::unique_ptr<IMarginalDistribution>> marginals;
                marginals.reserve(m_types.size());

                for (std::size_t i = 0; i < m_types.size(); ++i) {
                    unsigned int s = base_seed + static_cast<unsigned int>(i);
                    MarginalConfig cfg = to_marginal_config(m_cfgs[i]);
                    auto up = DistributionFactory::create(m_types[i], std::move(cfg), s);
                    marginals.push_back(std::move(up));
                }

                unsigned int cop_seed = base_seed + 1000003u; // to separate
                CopulaConfig ccfg = to_copula_config(c_cfg);
                auto copula = CopulaFactory::create(c_type, std::move(ccfg), cop_seed);

                return std::make_shared<JointDistribution>(std::move(marginals), std::move(copula));
            },
            py::arg("marginal_types"),
            py::arg("marginal_configs"),
            py::arg("copula_type"),
            py::arg("copula_config"),
            py::arg("seed") = py::none(),
            "Create a JointDistribution from (marginal types+configs) and (copula type+config)."
        )

        .def_static(
            "create_with_seeds",
            [](const std::vector<MarginalType> &m_types,
               const std::vector<py::object> &m_cfgs,
               const std::vector<unsigned int> &m_seeds,
               CopulaType c_type,
               py::object c_cfg,
               unsigned int copula_seed) {

                if (m_types.size() != m_cfgs.size() || m_types.size() != m_seeds.size()) {
                    throw py::value_error("marginal_types, marginal_configs, marginal_seeds must have the same length.");
                }
                if (m_types.empty()) {
                    throw py::value_error("At least one marginal is required.");
                }

                std::vector<std::unique_ptr<IMarginalDistribution>> marginals;
                marginals.reserve(m_types.size());

                for (std::size_t i = 0; i < m_types.size(); ++i) {
                    MarginalConfig cfg = to_marginal_config(m_cfgs[i]);
                    auto up = DistributionFactory::create(m_types[i], std::move(cfg), m_seeds[i]);
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

void init_statistic(py::module &m) {

    init_marginals(m);
    init_marginal_config_factory(m);
    init_copulas(m);
    init_joint_distribution(m);

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
    .def_readwrite("obss", &StatisticConfig::obss)
    .def_readwrite("p_specs", &StatisticConfig::p_specs)
    .def_readwrite("MC_draws", &StatisticConfig::MC_draws)
    .def_readwrite("skew_abs_threshold", &StatisticConfig::skew_abs_threshold)
    .def_readwrite("MLE_max_iter", &StatisticConfig::MLE_max_iter)
    .def_readwrite("MLE_tol", &StatisticConfig::MLE_tol);

    py::class_<FitResultWithMaps>(m, "FitResultWithMaps");

    py::class_<StatisticInterface, std::shared_ptr<StatisticInterface>>(m, "StatisticInterface")
     .def(py::init<StatisticConfig>(), py::arg("config"))

     .def("compute_uncertainties", &StatisticInterface::compute_uncertainties)
     .def("compute_uncertainties_and_sampling",
         &StatisticInterface::compute_uncertainties_and_sampling);
    //  .def("compute_MLE", &StatisticInterface::compute_MLE);
    
}

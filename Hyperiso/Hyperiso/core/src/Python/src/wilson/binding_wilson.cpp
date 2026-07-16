#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>

#include <utility>

#include "WilsonManager.h"
#include "MartyWilson.h"
#include "WilsonInterface.h"
#include "CustomWilsonLambda.h"
#include "SourcesView.h"
#include "BWilsonSUSY.h"
#include "BWilsonTHDM.h"

namespace py = pybind11;

namespace {

scalar_t scalar_from_python(py::handle obj) {
    if (py::hasattr(obj, "_cpp_obj")) {
        return py::getattr(obj, "_cpp_obj").cast<scalar_t>();
    }

    try {
        return obj.cast<scalar_t>();
    } catch (const py::cast_error&) {
    }

    try {
        return scalar_t(obj.cast<std::complex<double>>());
    } catch (const py::cast_error&) {
    }

    return scalar_t(obj.cast<double>());
}

std::function<scalar_t(const ParamSrc&)> wrap_matching_callable(py::function fn) {
    return [fn = std::move(fn)](const ParamSrc& src) -> scalar_t {
        py::gil_scoped_acquire gil;
        py::object out = fn(py::cast(src, py::return_value_policy::reference));
        return scalar_from_python(out);
    };
}

std::function<std::unordered_map<WCoefId, scalar_t>(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>&,
    const BlockSrc&
)> wrap_running_callable(py::function fn) {
    return [fn = std::move(fn)](
        const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& matching,
        const BlockSrc& src
    ) -> std::unordered_map<WCoefId, scalar_t> {
        py::gil_scoped_acquire gil;
        py::object out = fn(matching, py::cast(src, py::return_value_policy::reference));
        std::unordered_map<WCoefId, scalar_t> converted;
        py::dict d = out.cast<py::dict>();
        for (auto item : d) {
            converted.emplace(item.first.cast<WCoefId>(), scalar_from_python(item.second));
        }
        return converted;
    };
}

} // namespace

void init_coefficient_groups(py::module &m) {}

void init_wilson_coefficient(py::module &m) {}

void init_coefficient_manager(py::module &m) {}

void init_wilson_parameters(py::module &m) {}

void init_custom_wilson_lambda(py::module &m) {
    py::class_<CustomWilsonCoefficientConfig>(m, "CustomWilsonCoefficientConfig", R"pbdoc(
Runtime definition of a custom Wilson coefficient.

A coefficient config binds a dynamic ``WCoefId`` to one or more matching lambdas.
Each matching lambda declares the ``ParamId`` sources it needs, which makes the
coefficient usable by dependency-aware observable and statistic workflows.

The Python callback receives a ``ParamSrc`` view and must return a scalar value
(``float``/``complex`` depending on the scalar type used by Hyperiso).
)pbdoc")
        .def(py::init<>())
        .def(py::init<WCoefId>(), py::arg("id"))
        .def_readwrite("id", &CustomWilsonCoefficientConfig::id)
        .def("set_matching",
            [](CustomWilsonCoefficientConfig& self,
               QCDOrder order,
               std::unordered_set<ParamId> sources,
               py::function compute,
               ContributionType contribution) -> CustomWilsonCoefficientConfig& {
                return self.set_matching(
                    order,
                    std::move(sources),
                    wrap_matching_callable(std::move(compute)),
                    contribution
                );
            },
            py::arg("order"),
            py::arg("sources"),
            py::arg("compute"),
            py::arg("contribution") = ContributionType::SM,
            py::return_value_policy::reference_internal,
            R"pbdoc(Register a matching lambda for one QCD order.)pbdoc");

    py::class_<CustomWilsonGroupConfig>(m, "CustomWilsonGroupConfig", R"pbdoc(
Runtime definition of a Wilson group backed by user lambdas.

The group owns dynamic coefficients, matching/running scales and optional
running lambdas. If no running lambda is supplied, the C++ layer can install an
identity running in ``B_STANDARD`` so custom coefficients remain queryable at the
hadronic scale.
)pbdoc")
        .def(py::init<>())
        .def(py::init<WGroupId>(), py::arg("group"))
        .def_readwrite("group", &CustomWilsonGroupConfig::group)
        .def_readwrite("display_name", &CustomWilsonGroupConfig::display_name)
        .def_readwrite("matching_scale", &CustomWilsonGroupConfig::matching_scale)
        .def_readwrite("hadronic_scale", &CustomWilsonGroupConfig::hadronic_scale)
        .def_readwrite("order", &CustomWilsonGroupConfig::order)
        .def_readwrite("contribution", &CustomWilsonGroupConfig::contribution)
        .def_readwrite("coefficients", &CustomWilsonGroupConfig::coefficients)
        .def_readwrite("install_identity_running_if_empty", &CustomWilsonGroupConfig::install_identity_running_if_empty)
        .def("add_coefficient",
            [](CustomWilsonGroupConfig& self, CustomWilsonCoefficientConfig coef) -> CustomWilsonGroupConfig& {
                return self.add_coefficient(std::move(coef));
            },
            py::arg("coefficient"),
            py::return_value_policy::reference_internal,
            R"pbdoc(Append a coefficient config to this group.)pbdoc")
        .def("set_running",
            [](CustomWilsonGroupConfig& self,
               WilsonBasis basis,
               QCDOrder order,
               std::unordered_map<ParameterType, std::vector<std::string>> sources,
               py::function compute) -> CustomWilsonGroupConfig& {
                return self.set_running(
                    basis,
                    order,
                    std::move(sources),
                    wrap_running_callable(std::move(compute))
                );
            },
            py::arg("basis"),
            py::arg("order"),
            py::arg("sources"),
            py::arg("compute"),
            py::return_value_policy::reference_internal,
            R"pbdoc(Register a running lambda for one basis and QCD order.)pbdoc");
}

void init_wilson_interface(py::module &m) {
    init_custom_wilson_lambda(m);

    py::class_<WilsonInterface, std::shared_ptr<WilsonInterface>>(m, "WilsonInterface")
        .def(py::init<>())
        .def("build", &WilsonInterface::build, py::arg("wilson_config"), py::call_guard<py::gil_scoped_release>())
        .def("add_wilson_group", &WilsonInterface::addWilsonGroup, py::arg("config"), py::call_guard<py::gil_scoped_release>())
        .def("add_custom_group", &WilsonInterface::add_custom_group,
             py::arg("config"), py::return_value_policy::reference_internal,
             R"pbdoc(Add a lambda-backed custom Wilson group.)pbdoc")
        .def("addCustomWilsonGroup", &WilsonInterface::addCustomWilsonGroup,
             py::arg("config"), py::return_value_policy::reference_internal,
             R"pbdoc(CamelCase alias for add_custom_group.)pbdoc")
        .def("set_matching_scale", &WilsonInterface::set_matching_scale, py::arg("mu_W"))
        .def("set_hadronic_scale", &WilsonInterface::set_hadronic_scale, py::arg("mu_h"))
        .def("get_matching_coefficient", py::overload_cast<WGroup, WCoef, QCDOrder, ContributionType>(&WilsonInterface::getMatchingCoefficient))
        .def("get_matching_coefficient", py::overload_cast<WGroupId, WCoefId, QCDOrder, ContributionType>(&WilsonInterface::getMatchingCoefficient))
        .def("get_M", py::overload_cast<WGroup, WCoef, QCDOrder, ContributionType>(&WilsonInterface::getM))
        .def("get_M", py::overload_cast<WGroupId, WCoefId, QCDOrder, ContributionType>(&WilsonInterface::getM))
        .def("get_full_matching_coefficient", py::overload_cast<WGroup, WCoef, QCDOrder, ContributionType>(&WilsonInterface::getFullMatchingCoefficient))
        .def("get_full_matching_coefficient", py::overload_cast<WGroupId, WCoefId, QCDOrder, ContributionType>(&WilsonInterface::getFullMatchingCoefficient))
        .def("get_FM", py::overload_cast<WGroup, WCoef, QCDOrder, ContributionType>(&WilsonInterface::getFM))
        .def("get_FM", py::overload_cast<WGroupId, WCoefId, QCDOrder, ContributionType>(&WilsonInterface::getFM))
        .def("get_run_coefficient", py::overload_cast<WGroup, WCoef, QCDOrder, ContributionType, WilsonBasis>(&WilsonInterface::getRunCoefficient),
             py::arg("group"), py::arg("coeff"), py::arg("order"), py::arg("contribution"), py::arg("basis") = WilsonBasis::B_STANDARD)
        .def("get_run_coefficient", py::overload_cast<WGroupId, WCoefId, QCDOrder, ContributionType, WilsonBasis>(&WilsonInterface::getRunCoefficient),
             py::arg("group"), py::arg("coeff"), py::arg("order"), py::arg("contribution"), py::arg("basis") = WilsonBasis::B_STANDARD)
        .def("get_R", py::overload_cast<WGroup, WCoef, QCDOrder, ContributionType, WilsonBasis>(&WilsonInterface::getR),
             py::arg("group"), py::arg("coeff"), py::arg("order"), py::arg("contribution"), py::arg("basis") = WilsonBasis::B_STANDARD)
        .def("get_R", py::overload_cast<WGroupId, WCoefId, QCDOrder, ContributionType, WilsonBasis>(&WilsonInterface::getR),
             py::arg("group"), py::arg("coeff"), py::arg("order"), py::arg("contribution"), py::arg("basis") = WilsonBasis::B_STANDARD)
        .def("get_full_run_coefficient", py::overload_cast<WGroup, WCoef, QCDOrder, ContributionType, WilsonBasis>(&WilsonInterface::getFullRunCoefficient),
             py::arg("group"), py::arg("coeff"), py::arg("order"), py::arg("contribution"), py::arg("basis") = WilsonBasis::B_STANDARD)
        .def("get_full_run_coefficient", py::overload_cast<WGroupId, WCoefId, QCDOrder, ContributionType, WilsonBasis>(&WilsonInterface::getFullRunCoefficient),
             py::arg("group"), py::arg("coeff"), py::arg("order"), py::arg("contribution"), py::arg("basis") = WilsonBasis::B_STANDARD)
        .def("get_FR", py::overload_cast<WGroup, WCoef, QCDOrder, ContributionType, WilsonBasis>(&WilsonInterface::getFR),
             py::arg("group"), py::arg("coeff"), py::arg("order"), py::arg("contribution"), py::arg("basis") = WilsonBasis::B_STANDARD)
        .def("get_FR", py::overload_cast<WGroupId, WCoefId, QCDOrder, ContributionType, WilsonBasis>(&WilsonInterface::getFR),
             py::arg("group"), py::arg("coeff"), py::arg("order"), py::arg("contribution"), py::arg("basis") = WilsonBasis::B_STANDARD)
        .def("get_sep_order_matching_coefficient", py::overload_cast<WGroup, WCoef, ContributionType>(&WilsonInterface::getSepOrderMatchingCoefficient))
        .def("get_sep_order_matching_coefficient", py::overload_cast<WGroupId, WCoefId, ContributionType>(&WilsonInterface::getSepOrderMatchingCoefficient))
        .def("get_sep_order_run_coefficient", py::overload_cast<WGroup, WCoef, ContributionType, WilsonBasis>(&WilsonInterface::getSepOrderRunCoefficient),
             py::arg("group"), py::arg("coeff"), py::arg("contribution"), py::arg("basis") = WilsonBasis::B_STANDARD)
        .def("get_sep_order_run_coefficient", py::overload_cast<WGroupId, WCoefId, ContributionType, WilsonBasis>(&WilsonInterface::getSepOrderRunCoefficient),
             py::arg("group"), py::arg("coeff"), py::arg("contribution"), py::arg("basis") = WilsonBasis::B_STANDARD)
        .def("get_SM", py::overload_cast<WGroup, WCoef, ContributionType>(&WilsonInterface::getSM))
        .def("get_SM", py::overload_cast<WGroupId, WCoefId, ContributionType>(&WilsonInterface::getSM))
        .def("get_SR", py::overload_cast<WGroup, WCoef, ContributionType, WilsonBasis>(&WilsonInterface::getSR),
             py::arg("group"), py::arg("coeff"), py::arg("contribution"), py::arg("basis") = WilsonBasis::B_STANDARD)
        .def("get_SR", py::overload_cast<WGroupId, WCoefId, ContributionType, WilsonBasis>(&WilsonInterface::getSR),
             py::arg("group"), py::arg("coeff"), py::arg("contribution"), py::arg("basis") = WilsonBasis::B_STANDARD)
        .def("get_all_matching_coefficient", &WilsonInterface::getAllMatchingCoefficients)
        .def("get_all_run_coefficient", &WilsonInterface::getAllRunCoefficients)
        .def("get_all_full_matching_coefficient", &WilsonInterface::getAllFullMatchingCoefficients)
        .def("get_all_full_run_coefficient", &WilsonInterface::getAllFullRunCoefficients);
}

void init_wilson(py::module &m) {
    auto coeff_groups = m.def_submodule("coefficient_groups", "Coefficient group management");
    init_coefficient_groups(coeff_groups);

    auto coeff = m.def_submodule("wilson_coefficient", "Wilson coefficient management");
    init_wilson_coefficient(coeff);

    auto manager = m.def_submodule("coefficient_manager", "Coefficient manager");
    init_coefficient_manager(manager);

    auto params = m.def_submodule("wilson_parameters", "Wilson parameters");
    init_wilson_parameters(params);

    auto interface = m.def_submodule("wilson_interface", "Wilson interface management");
    init_wilson_interface(interface);
}

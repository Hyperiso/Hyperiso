#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>

#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>

#include "ObservableInterface.h"


#include "DefaultConfig.h"
#include "BDlnuDecay.h"
#include "BDstarlnuDecay.h"
#include "BKllDecay.h"
#include "BKsllDecay.h"
#include "BKstarGammaDecay.h"
#include "BllDecay.h"
#include "BsPhiDecay.h"
#include "BXsDecay.h"
#include "BXsllDecay.h"
#include "DlnuDecay.h"
#include "DslnuDecay.h"
#include "KllDecay.h"
#include "KlnuDecay.h"
#include "KPinunuDecay.h"
#include "LbLllDecay.h"
#include "M0_Mixing.h"

namespace py = pybind11;

namespace {

double double_from_python(py::handle obj) {
    if (py::hasattr(obj, "_cpp_obj")) {
        py::object inner = py::getattr(obj, "_cpp_obj");
        try {
            return inner.cast<double>();
        } catch (const py::cast_error&) {
            // Continue with complex-like values.
        }
        try {
            return inner.cast<std::complex<double>>().real();
        } catch (const py::cast_error&) {
            // Continue with the public object itself.
        }
    }

    try {
        return obj.cast<double>();
    } catch (const py::cast_error&) {
        // Continue with Python complex.
    }

    return obj.cast<std::complex<double>>().real();
}

std::function<double(LambdaDecay&, ObservableId)> wrap_lambda_scalar(py::function fn) {
    return [fn = std::move(fn)](LambdaDecay& ctx, ObservableId id) -> double {
        py::gil_scoped_acquire gil;
        py::object out = fn(py::cast(&ctx, py::return_value_policy::reference), id);
        return double_from_python(out);
    };
}

std::function<double(LambdaDecay&, std::pair<double,double>, ObservableId)> wrap_lambda_binned_scalar(py::function fn) {
    return [fn = std::move(fn)](LambdaDecay& ctx, std::pair<double,double> bin, ObservableId id) -> double {
        py::gil_scoped_acquire gil;
        py::object out = fn(py::cast(&ctx, py::return_value_policy::reference), bin, id);
        return double_from_python(out);
    };
}

void bind_lambda_decay_types(py::module& m) {
    py::class_<LambdaDecay>(m, "LambdaDecay", R"pbdoc(
Execution context passed to Python lambda observables.

The context exposes the same services that a compiled ``DecayParent`` uses:
Wilson coefficients, SM/FLAVOR parameters, QCD helpers and the current bin list.
The object is owned by the C++ observable manager; keep it only during the
callback call.
)pbdoc")
        .def("get_M", [](LambdaDecay& ctx, WGroupId group, WCoefId coeff, QCDOrder order, ContributionType contribution) {
            return ctx.W().getM(group, coeff, order, contribution);
        }, py::arg("group"), py::arg("coeff"), py::arg("order"), py::arg("contribution"))
        .def("get_FM", [](LambdaDecay& ctx, WGroupId group, WCoefId coeff, QCDOrder order, ContributionType contribution) {
            return ctx.W().getFM(group, coeff, order, contribution);
        }, py::arg("group"), py::arg("coeff"), py::arg("order"), py::arg("contribution"))
        .def("get_R", [](LambdaDecay& ctx, WGroupId group, WCoefId coeff, QCDOrder order, ContributionType contribution) {
            return ctx.W().getR(group, coeff, order, contribution);
        }, py::arg("group"), py::arg("coeff"), py::arg("order"), py::arg("contribution"))
        .def("get_FR", [](LambdaDecay& ctx, WGroupId group, WCoefId coeff, QCDOrder order, ContributionType contribution) {
            return ctx.W().getFR(group, coeff, order, contribution);
        }, py::arg("group"), py::arg("coeff"), py::arg("order"), py::arg("contribution"))
        .def("get_sm_param", [](LambdaDecay& ctx, const ParamId& pid, DataType type) {
            return ctx.SM()(pid, type);
        }, py::arg("pid"), py::arg("type") = DataType::VALUE)
        .def("get_flavor_param", [](LambdaDecay& ctx, const ParamId& pid, DataType type) {
            return ctx.FLAVOR()(pid, type);
        }, py::arg("pid"), py::arg("type") = DataType::VALUE)
        .def("current_bins", [](LambdaDecay& ctx) {
            auto bins = ctx.current_bins();
            if (!bins.has_value()) {
                return std::vector<std::pair<double,double>>{};
            }
            return *bins;
        }, R"pbdoc(Return currently requested bins, or an empty list for unbinned calls.)pbdoc");

    py::class_<LambdaObservableConfig>(m, "LambdaObservableConfig", R"pbdoc(
Runtime observable definition for ``LambdaDecayConfig``.

Use ``scalar`` for unbinned observables and ``binned_scalar`` for observables
that should be evaluated once per currently selected bin. Dependencies declared
on this object are visible to ``ObservableInterface`` and ``StatisticInterface``.
)pbdoc")
        .def(py::init<>())
        .def_readwrite("canonical", &LambdaObservableConfig::canonical)
        .def_readwrite("aliases", &LambdaObservableConfig::aliases)
        .def_readwrite("flha", &LambdaObservableConfig::flha)
        .def_readwrite("dependencies", &LambdaObservableConfig::dependencies)
        .def_static("scalar", [](std::string canonical, py::function compute) {
            return LambdaObservableConfig::scalar(std::move(canonical), wrap_lambda_scalar(std::move(compute)));
        }, py::arg("canonical"), py::arg("compute"),
           R"pbdoc(Create an unbinned scalar observable from a Python callback.)pbdoc")
        .def_static("binned_scalar", [](std::string canonical, py::function compute) {
            return LambdaObservableConfig::binned_scalar(std::move(canonical), wrap_lambda_binned_scalar(std::move(compute)));
        }, py::arg("canonical"), py::arg("compute"),
           R"pbdoc(Create a binned scalar observable from a Python callback.)pbdoc");

    py::class_<LambdaDecayConfig>(m, "LambdaDecayConfig", R"pbdoc(
Runtime decay definition backed by Python lambdas.

The config registers a new dynamic decay, optional custom Wilson groups and one
or more lambda observables. When ``propagate_custom_wilson_dependencies`` is
true, ParamId sources declared on custom Wilson matching lambdas are attached to
the observables so the statistic layer can discover them.
)pbdoc")
        .def(py::init<>())
        .def_readwrite("canonical", &LambdaDecayConfig::canonical)
        .def_readwrite("aliases", &LambdaDecayConfig::aliases)
        .def_readwrite("matching_scale", &LambdaDecayConfig::matching_scale)
        .def_readwrite("hadronic_scale", &LambdaDecayConfig::hadronic_scale)
        .def_readwrite("order", &LambdaDecayConfig::order)
        .def_readwrite("max_order", &LambdaDecayConfig::max_order)
        .def_readwrite("wilson_groups", &LambdaDecayConfig::wilson_groups)
        .def_readwrite("custom_wilson_groups", &LambdaDecayConfig::custom_wilson_groups)
        .def_readwrite("observables", &LambdaDecayConfig::observables)
        .def_readwrite("propagate_custom_wilson_dependencies", &LambdaDecayConfig::propagate_custom_wilson_dependencies);
}

void bind_decay_config_types(py::module& m) {
    py::class_<DecayConfig>(m, "DecayConfig", R"pbdoc(
Base decay-configuration marker used by configurable decay implementations.

This type has no fields by itself. Some concrete decay configuration objects
inherit from it in C++, while other legacy config structs are stored directly in
``std::any``. All of them are passed to ``ObservableInterface.set_decay_config``
through typed pybind11 overloads.
)pbdoc")
        .def(py::init<>(), R"pbdoc(Create an empty base decay configuration.)pbdoc");

    py::enum_<B_FF_Type>(m, "B_FF_Type", R"pbdoc(
Choice of form-factor scheme for exclusive B decays.
)pbdoc")
        .value("FULL", B_FF_Type::FULL, R"pbdoc(Use full form factors.)pbdoc")
        .value("SOFT", B_FF_Type::SOFT, R"pbdoc(Use soft form-factor approximation.)pbdoc");

    py::enum_<BP_FF_Src>(m, "BP_FF_Src", R"pbdoc(
Source of B -> pseudoscalar form factors used by BKllConfig.
)pbdoc")
        .value("AS", BP_FF_Src::AS)
        .value("GRvDV", BP_FF_Src::GRvDV)
        .value("GKvD_SR_LAT", BP_FF_Src::GKvD_SR_LAT)
        .value("GKvD_SR", BP_FF_Src::GKvD_SR)
        .value("FLAG24", BP_FF_Src::FLAG24)
        .value("HPQCD22", BP_FF_Src::HPQCD22);

    py::enum_<BV_FF_Src>(m, "BV_FF_Src", R"pbdoc(
Source of B -> vector form factors used by BKstarllConfig, BKstarGammaConfig and BsPhiConfig.
)pbdoc")
        .value("BSZ_SR_LAT", BV_FF_Src::BSZ_SR_LAT)
        .value("BSZ_SR", BV_FF_Src::BSZ_SR)
        .value("GRvDV", BV_FF_Src::GRvDV)
        .value("GKvD_SR_LAT", BV_FF_Src::GKvD_SR_LAT)
        .value("GKvD_SR", BV_FF_Src::GKvD_SR)
        .value("HLMW", BV_FF_Src::HLMW);

    py::enum_<LbL_FF_Src>(m, "LbL_FF_Src", R"pbdoc(
Source of Lambda_b -> Lambda form factors used by LbLllConfig.

Only values available in the current C++ headers are exposed here.
)pbdoc")
        .value("DM", LbL_FF_Src::DM);

    auto bdlnu_cfg = py::class_<BDlnuConfig, DecayConfig>(m, "BDlnuConfig", R"pbdoc(
Configuration for the B -> D l nu decay engine.

Attributes:
    charge: B-meson charge convention used by the decay implementation.
)pbdoc");
    py::enum_<BDlnuConfig::B_Charge>(bdlnu_cfg, "BCharge", R"pbdoc(B charge option for BDlnuConfig.)pbdoc")
        .value("B_0", BDlnuConfig::B_Charge::B_0)
        .value("B_PLUS", BDlnuConfig::B_Charge::B_PLUS);
    bdlnu_cfg
        .def(py::init<>(), R"pbdoc(Create a B -> D l nu configuration with backend defaults.)pbdoc")
        .def_readwrite("charge", &BDlnuConfig::charge, R"pbdoc(B-meson charge option.)pbdoc");

    auto bdstarlnu_cfg = py::class_<BDstarlnuConfig, DecayConfig>(m, "BDstarlnuConfig", R"pbdoc(
Configuration for the B -> D* l nu decay engine.

Attributes:
    charge: B-meson charge convention used by the decay implementation.
)pbdoc");
    py::enum_<BDstarlnuConfig::B_Charge>(bdstarlnu_cfg, "BCharge", R"pbdoc(B charge option for BDstarlnuConfig.)pbdoc")
        .value("B_0", BDstarlnuConfig::B_Charge::B_0)
        .value("B_PLUS", BDstarlnuConfig::B_Charge::B_PLUS);
    bdstarlnu_cfg
        .def(py::init<>(), R"pbdoc(Create a B -> D* l nu configuration with backend defaults.)pbdoc")
        .def_readwrite("charge", &BDstarlnuConfig::charge, R"pbdoc(B-meson charge option.)pbdoc");

    auto bkll_cfg = py::class_<BKllConfig, DecayConfig>(m, "BKllConfig", R"pbdoc(
Configuration for the exclusive B -> K l+ l- decay engine.

Attributes:
    ff_src: Source of B -> K form factors.
    ff_type: Full or soft form-factor treatment.
    charge: B-meson charge convention.
    gen: Lepton generation.
    n_threads: Number of worker threads requested by the decay implementation.
)pbdoc");
    py::enum_<BKllConfig::B_Charge>(bkll_cfg, "BCharge", R"pbdoc(B charge option for BKllConfig.)pbdoc")
        .value("B_0", BKllConfig::B_Charge::B_0)
        .value("B_PLUS", BKllConfig::B_Charge::B_PLUS);
    py::enum_<BKllConfig::Lepton>(bkll_cfg, "Lepton", R"pbdoc(Lepton generation option for BKllConfig.)pbdoc")
        .value("E", BKllConfig::Lepton::E)
        .value("MU", BKllConfig::Lepton::MU)
        .value("TAU", BKllConfig::Lepton::TAU);
    bkll_cfg
        .def(py::init<>(), R"pbdoc(Create a B -> K l+ l- configuration with backend defaults.)pbdoc")
        .def_readwrite("ff_src", &BKllConfig::ff_src, R"pbdoc(B -> pseudoscalar form-factor source.)pbdoc")
        .def_readwrite("ff_type", &BKllConfig::ff_type, R"pbdoc(Form-factor treatment, full or soft.)pbdoc")
        .def_readwrite("charge", &BKllConfig::charge, R"pbdoc(B-meson charge option.)pbdoc")
        .def_readwrite("gen", &BKllConfig::gen, R"pbdoc(Lepton generation.)pbdoc")
        .def_readwrite("n_threads", &BKllConfig::n_threads, R"pbdoc(Number of worker threads.)pbdoc");

    auto bkstarll_cfg = py::class_<BKstarllConfig, DecayConfig>(m, "BKstarllConfig", R"pbdoc(
Configuration for the exclusive B -> K* l+ l- decay engine.

Attributes:
    ff_src: Source of B -> K* form factors.
    ff_type: Full or soft form-factor treatment.
    power_corr_impl: Non-factorisable power-correction prescription.
    charge: B-meson charge convention.
    gen: Lepton generation.
    n_threads: Number of worker threads requested by the decay implementation.
)pbdoc");
    py::enum_<BKstarllConfig::Power_Corrections_Impl>(bkstarll_cfg, "PowerCorrectionsImpl", R"pbdoc(Power-correction model for BKstarllConfig.)pbdoc")
        .value("BFS", BKstarllConfig::Power_Corrections_Impl::BFS)
        .value("BCvDV", BKstarllConfig::Power_Corrections_Impl::BCvDV)
        .value("KMPW", BKstarllConfig::Power_Corrections_Impl::KMPW);
    py::enum_<BKstarllConfig::B_Charge>(bkstarll_cfg, "BCharge", R"pbdoc(B charge option for BKstarllConfig.)pbdoc")
        .value("B_0", BKstarllConfig::B_Charge::B_0)
        .value("B_PLUS", BKstarllConfig::B_Charge::B_PLUS);
    py::enum_<BKstarllConfig::Lepton>(bkstarll_cfg, "Lepton", R"pbdoc(Lepton generation option for BKstarllConfig.)pbdoc")
        .value("E", BKstarllConfig::Lepton::E)
        .value("MU", BKstarllConfig::Lepton::MU)
        .value("TAU", BKstarllConfig::Lepton::TAU);
    bkstarll_cfg
        .def(py::init<>(), R"pbdoc(Create a B -> K* l+ l- configuration with backend defaults.)pbdoc")
        .def_readwrite("ff_src", &BKstarllConfig::ff_src, R"pbdoc(B -> vector form-factor source.)pbdoc")
        .def_readwrite("ff_type", &BKstarllConfig::ff_type, R"pbdoc(Form-factor treatment, full or soft.)pbdoc")
        .def_readwrite("power_corr_impl", &BKstarllConfig::power_corr_impl, R"pbdoc(Power-correction prescription.)pbdoc")
        .def_readwrite("charge", &BKstarllConfig::charge, R"pbdoc(B-meson charge option.)pbdoc")
        .def_readwrite("gen", &BKstarllConfig::gen, R"pbdoc(Lepton generation.)pbdoc")
        .def_readwrite("n_threads", &BKstarllConfig::n_threads, R"pbdoc(Number of worker threads.)pbdoc");

    auto bkstargamma_cfg = py::class_<BKstarGammaConfig, DecayConfig>(m, "BKstarGammaConfig", R"pbdoc(
Configuration for the exclusive B -> K* gamma decay engine.

Attributes:
    ff_src: Source of B -> K* form factors.
    charge: B-meson charge convention.
)pbdoc");
    py::enum_<BKstarGammaConfig::B_Charge>(bkstargamma_cfg, "BCharge", R"pbdoc(B charge option for BKstarGammaConfig.)pbdoc")
        .value("B_0", BKstarGammaConfig::B_Charge::B_0)
        .value("B_PLUS", BKstarGammaConfig::B_Charge::B_PLUS);
    bkstargamma_cfg
        .def(py::init<>(), R"pbdoc(Create a B -> K* gamma configuration with backend defaults.)pbdoc")
        .def_readwrite("ff_src", &BKstarGammaConfig::ff_src, R"pbdoc(B -> vector form-factor source.)pbdoc")
        .def_readwrite("charge", &BKstarGammaConfig::charge, R"pbdoc(B-meson charge option.)pbdoc");

    auto bsphi_cfg = py::class_<BsPhiConfig>(m, "BsPhiConfig", R"pbdoc(
Configuration for the exclusive Bs -> phi l+ l- decay engine.

Attributes:
    ff_src: Source of Bs -> phi form factors.
    ff_type: Full or soft form-factor treatment.
    gen: Lepton generation.
    n_threads: Number of worker threads requested by the decay implementation.
)pbdoc");
    py::enum_<BsPhiConfig::Lepton>(bsphi_cfg, "Lepton", R"pbdoc(Lepton generation option for BsPhiConfig.)pbdoc")
        .value("E", BsPhiConfig::Lepton::E)
        .value("MU", BsPhiConfig::Lepton::MU)
        .value("TAU", BsPhiConfig::Lepton::TAU);
    bsphi_cfg
        .def(py::init<>(), R"pbdoc(Create a Bs -> phi l+ l- configuration with backend defaults.)pbdoc")
        .def_readwrite("ff_src", &BsPhiConfig::ff_src, R"pbdoc(Bs -> phi form-factor source.)pbdoc")
        .def_readwrite("ff_type", &BsPhiConfig::ff_type, R"pbdoc(Form-factor treatment, full or soft.)pbdoc")
        .def_readwrite("gen", &BsPhiConfig::gen, R"pbdoc(Lepton generation.)pbdoc")
        .def_readwrite("n_threads", &BsPhiConfig::n_threads, R"pbdoc(Number of worker threads.)pbdoc");

    auto bxsll_cfg = py::class_<BXsllConfig, DecayConfig>(m, "BXsllConfig", R"pbdoc(
Configuration for the inclusive B -> X_s l+ l- decay engine.

Attributes:
    gen: Lepton generation.
)pbdoc");
    py::enum_<BXsllConfig::Lepton>(bxsll_cfg, "Lepton", R"pbdoc(Lepton generation option for BXsllConfig.)pbdoc")
        .value("E", BXsllConfig::Lepton::E)
        .value("MU", BXsllConfig::Lepton::MU)
        .value("TAU", BXsllConfig::Lepton::TAU);
    bxsll_cfg
        .def(py::init<>(), R"pbdoc(Create a B -> X_s l+ l- configuration with backend defaults.)pbdoc")
        .def_readwrite("gen", &BXsllConfig::gen, R"pbdoc(Lepton generation.)pbdoc");

    py::class_<KllDecayConfig>(m, "KllDecayConfig", R"pbdoc(
Configuration for the K_L,S -> l+ l- decay engine.

Attributes:
    N_L_sign: Sign convention for the K_L -> gamma gamma long-distance term.
    gen: Lepton generation encoded as the integer convention used by the C++ backend.
)pbdoc")
        .def(py::init<>(), R"pbdoc(Create a K -> l+ l- configuration with backend defaults.)pbdoc")
        .def_readwrite("N_L_sign", &KllDecayConfig::N_L_sign, R"pbdoc(Sign convention for the long-distance K_L contribution.)pbdoc")
        .def_readwrite("gen", &KllDecayConfig::gen, R"pbdoc(Lepton-generation integer used by the backend.)pbdoc");

    auto lblll_cfg = py::class_<LbLllConfig>(m, "LbLllConfig", R"pbdoc(
Configuration for the Lambda_b -> Lambda l+ l- decay engine.

Attributes:
    ff_src: Source of Lambda_b -> Lambda form factors.
    gen: Lepton generation.
)pbdoc");
    py::enum_<LbLllConfig::Lepton>(lblll_cfg, "Lepton", R"pbdoc(Lepton generation option for LbLllConfig.)pbdoc")
        .value("E", LbLllConfig::Lepton::E)
        .value("MU", LbLllConfig::Lepton::MU)
        .value("TAU", LbLllConfig::Lepton::TAU);
    lblll_cfg
        .def(py::init<>(), R"pbdoc(Create a Lambda_b -> Lambda l+ l- configuration with backend defaults.)pbdoc")
        .def_readwrite("ff_src", &LbLllConfig::ff_src, R"pbdoc(Lambda_b -> Lambda form-factor source.)pbdoc")
        .def_readwrite("gen", &LbLllConfig::gen, R"pbdoc(Lepton generation.)pbdoc");
}

template <typename ConfigT>
ObservableInterface& set_decay_config_typed(ObservableInterface& self, Decays decay, const ConfigT& config) {
    self.set_decay_config(decay, config);
    return self;
}

} // namespace




void init_observable(py::module &m) {

    bind_decay_config_types(m);
    bind_lambda_decay_types(m);

    py::class_<ObservableValue>(m, "ObservableValue")
        .def(py::init<ObservableId, double>(),
             py::arg("id"), py::arg("value"))
        .def(py::init<ObservableId, double, std::pair<double,double>>(),
             py::arg("id"), py::arg("value"), py::arg("bin"))

        .def_readwrite("id", &ObservableValue::id)
        .def_readwrite("value", &ObservableValue::value)

        .def_property(
            "bin",
            [](const ObservableValue &ov) -> py::object {
                if (ov.bin) return py::cast(*ov.bin);
                return py::none();
            },
            [](ObservableValue &ov, py::object b) {
                if (b.is_none()) ov.bin.reset();
                else ov.bin = b.cast<std::pair<double,double>>();
            })

        .def("__repr__", [](const ObservableValue &ov) {
            std::ostringstream oss;
            oss << "ObservableValue(id=";
            oss << py::str(py::cast(ov.id)).cast<std::string>();
            oss << ", value=" << ov.value;
            if (ov.bin) oss << ", bin=(" << ov.bin->first << ", " << ov.bin->second << ")";
            oss << ")";
            return oss.str();
        });

    py::class_<ObservableInterface, std::shared_ptr<ObservableInterface>>(m, "ObservableInterface")
     .def(py::init<>())
     .def("add_lambda_decay", &ObservableInterface::add_lambda_decay,
          py::arg("config"), py::arg("add_observables") = true,
          py::return_value_policy::reference_internal,
          R"pbdoc(Register a lambda-backed custom decay and optionally add its observables.)pbdoc")

     // add_observable( … )
     .def("add_observable",
          py::overload_cast<Observables, QCDOrder, bool>(&ObservableInterface::add_observable),
          py::arg("obs"), py::arg("order"), py::arg("add_dependencies") = false)
     .def("add_observable",
          py::overload_cast<ObservableId, QCDOrder, bool>(&ObservableInterface::add_observable),
          py::arg("obs"), py::arg("order"), py::arg("add_dependencies") = false)
     .def("add_observable",
          py::overload_cast<BinnedObservableId, QCDOrder, bool>(&ObservableInterface::add_observable),
          py::arg("obs"), py::arg("order"), py::arg("add_dependencies") = false)

     // add_observables( … ) – exposes all public overloads except custom-decay/std::any APIs
     .def("add_observables", 
          py::overload_cast<std::map<Observables, QCDOrder>, bool>(&ObservableInterface::add_observables),
          py::arg("obss"), py::arg("add_dependencies") = false)
     .def("add_observables", 
          py::overload_cast<std::map<ObservableId, QCDOrder>, bool>(&ObservableInterface::add_observables),
          py::arg("obss"), py::arg("add_dependencies") = false)
     .def("add_observables", 
          py::overload_cast<Decays, QCDOrder, bool>(&ObservableInterface::add_observables),
          py::arg("decay"), py::arg("order"), py::arg("add_dependencies") = false)

     .def("add_observable_parameter",
          py::overload_cast<Observables, ParamId>(&ObservableInterface::add_observable_parameter),
          py::arg("obs"), py::arg("pid"))
     .def("add_observable_parameter",
          py::overload_cast<ObservableId, ParamId>(&ObservableInterface::add_observable_parameter),
          py::arg("obs"), py::arg("pid"))

     .def("add_observable_parameters",
          py::overload_cast<Observables, std::unordered_set<ParamId>>(&ObservableInterface::add_observable_parameters),
          py::arg("obs"), py::arg("pids"))
     .def("add_observable_parameters",
          py::overload_cast<ObservableId, std::unordered_set<ParamId>>(&ObservableInterface::add_observable_parameters),
          py::arg("obs"), py::arg("pids"))

     // compute_observable(…) const
     .def("compute_observable",
          py::overload_cast<Observables>(&ObservableInterface::compute_observable, py::const_),
          py::arg("obs"))
     .def("compute_observable",
          py::overload_cast<ObservableId>(&ObservableInterface::compute_observable, py::const_),
          py::arg("obs"))
     .def("compute_all", &ObservableInterface::compute_all)

     .def("remove_observable",
          py::overload_cast<Observables>(&ObservableInterface::remove_observable),
          py::arg("obs"))
     .def("remove_observable",
          py::overload_cast<ObservableId>(&ObservableInterface::remove_observable),
          py::arg("obs"))

     .def("remove_observables",
          py::overload_cast<std::unordered_set<Observables>>(
               &ObservableInterface::remove_observables),
          py::arg("ids"))
     .def("remove_observables",
          py::overload_cast<Decays>(
               &ObservableInterface::remove_observables),
          py::arg("decay"))

     // get_exp_value / get_exp_uncertainty
     .def("get_exp_value",
          py::overload_cast<Observables>(&ObservableInterface::get_exp_value),
          py::arg("obs"))
     .def("get_exp_value",
          py::overload_cast<ObservableId>(&ObservableInterface::get_exp_value),
          py::arg("obs"))

     .def("get_exp_uncertainty",
          py::overload_cast<Observables, UncertaintyType>(&ObservableInterface::get_exp_uncertainty),
          py::arg("obs"), py::arg("u_type") = UncertaintyType::COMBINED)
     .def("get_exp_uncertainty",
          py::overload_cast<ObservableId, UncertaintyType>(&ObservableInterface::get_exp_uncertainty),
          py::arg("obs"), py::arg("u_type") = UncertaintyType::COMBINED)

     .def("get_current_observables", &ObservableInterface::get_current_observables)
     .def("get_all_ops_deps",
          py::overload_cast<ObservableId>(&ObservableInterface::get_all_ops_deps),
          py::arg("obs"))
     .def("get_all_ops_deps",
          py::overload_cast<Observables>(&ObservableInterface::get_all_ops_deps),
          py::arg("obs"))
     .def("set_param", &ObservableInterface::set_param,
          py::arg("block"), py::arg("code"), py::arg("value"), py::arg("type"))
     .def("get_param", &ObservableInterface::get_param,
          py::arg("block"), py::arg("code"), py::arg("type"))
     .def("reload_params", &ObservableInterface::reload_params)
     .def("enable_obs", &ObservableInterface::enable_obs)

     .def("set_decay_config",
          [](ObservableInterface& self, Decays decay, const BDlnuConfig& config) -> ObservableInterface& {
              return set_decay_config_typed(self, decay, config);
          },
          py::arg("decay"), py::arg("config"), py::return_value_policy::reference_internal,
          R"pbdoc(Set the B -> D l nu decay configuration.)pbdoc")
     .def("set_decay_config",
          [](ObservableInterface& self, Decays decay, const BDstarlnuConfig& config) -> ObservableInterface& {
              return set_decay_config_typed(self, decay, config);
          },
          py::arg("decay"), py::arg("config"), py::return_value_policy::reference_internal,
          R"pbdoc(Set the B -> D* l nu decay configuration.)pbdoc")
     .def("set_decay_config",
          [](ObservableInterface& self, Decays decay, const BKllConfig& config) -> ObservableInterface& {
              return set_decay_config_typed(self, decay, config);
          },
          py::arg("decay"), py::arg("config"), py::return_value_policy::reference_internal,
          R"pbdoc(Set the B -> K l+ l- decay configuration.)pbdoc")
     .def("set_decay_config",
          [](ObservableInterface& self, Decays decay, const BKstarllConfig& config) -> ObservableInterface& {
              return set_decay_config_typed(self, decay, config);
          },
          py::arg("decay"), py::arg("config"), py::return_value_policy::reference_internal,
          R"pbdoc(Set the B -> K* l+ l- decay configuration.)pbdoc")
     .def("set_decay_config",
          [](ObservableInterface& self, Decays decay, const BKstarGammaConfig& config) -> ObservableInterface& {
              return set_decay_config_typed(self, decay, config);
          },
          py::arg("decay"), py::arg("config"), py::return_value_policy::reference_internal,
          R"pbdoc(Set the B -> K* gamma decay configuration.)pbdoc")
     .def("set_decay_config",
          [](ObservableInterface& self, Decays decay, const BsPhiConfig& config) -> ObservableInterface& {
              return set_decay_config_typed(self, decay, config);
          },
          py::arg("decay"), py::arg("config"), py::return_value_policy::reference_internal,
          R"pbdoc(Set the Bs -> phi l+ l- decay configuration.)pbdoc")
     .def("set_decay_config",
          [](ObservableInterface& self, Decays decay, const BXsllConfig& config) -> ObservableInterface& {
              return set_decay_config_typed(self, decay, config);
          },
          py::arg("decay"), py::arg("config"), py::return_value_policy::reference_internal,
          R"pbdoc(Set the B -> X_s l+ l- decay configuration.)pbdoc")
     .def("set_decay_config",
          [](ObservableInterface& self, Decays decay, const KllDecayConfig& config) -> ObservableInterface& {
              return set_decay_config_typed(self, decay, config);
          },
          py::arg("decay"), py::arg("config"), py::return_value_policy::reference_internal,
          R"pbdoc(Set the K -> l+ l- decay configuration.)pbdoc")
     .def("set_decay_config",
          [](ObservableInterface& self, Decays decay, const LbLllConfig& config) -> ObservableInterface& {
              return set_decay_config_typed(self, decay, config);
          },
          py::arg("decay"), py::arg("config"), py::return_value_policy::reference_internal,
          R"pbdoc(Set the Lambda_b -> Lambda l+ l- decay configuration.)pbdoc")
     .def("set_decay_config",
          [](ObservableInterface& self, Decays decay, const DecayConfig& config) -> ObservableInterface& {
              return set_decay_config_typed(self, decay, config);
          },
          py::arg("decay"), py::arg("config"), py::return_value_policy::reference_internal,
          R"pbdoc(Set a generic, empty decay configuration for a configurable decay.)pbdoc")
     .def("set_bkstarll_threads", &ObservableInterface::set_bkstarll_threads,
          py::arg("n_threads"))
     .def("set_bkll_threads", &ObservableInterface::set_bkll_threads,
          py::arg("n_threads"))
     .def("set_bsphi_threads", &ObservableInterface::set_bsphi_threads,
          py::arg("n_threads"));

}
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <optional>
#include <filesystem>
#include "Include.h"

#include "HyperisoMaster.h"
#include "Config.h"
#include "ParameterSetter.h"
#include "ParameterProvider.h"
#include "CorrelationProvider.h"
#include "QCDProvider.h"
#include "Freezer.h"
#include "DBMemento.h"
#include "APIAdapter.h"
#include "ParameterShifter.h"
#include "MartyAdapter.h"
#include "BlockProvider.h"
#include "QEDProvider.h"
#include "DependantBlockInfoProvider.h"

namespace py = pybind11;

template <typename T>
void declare_parameter(py::module &m, const std::string &name) {
    py::class_<T, std::shared_ptr<T>>(m, name.c_str())
        .def(py::init<>())
        .def(py::init<ParamId, scalar_t, scalar_t, scalar_t>(), py::arg("id"), py::arg("mean"), py::arg("std_stat"), py::arg("std_syst"))
        .def("get_val", &T::get_val)
        .def("get_combined_std", &T::get_combined_std)
        .def("get_std", &T::get_std)
        .def("get_id", &T::get_id)
        .def("get_scale", &T::get_scale)
        .def("get_bin", &T::get_bin)
        .def("set_expected", &T::set_expected)
        .def("set_std", &T::set_std)
        .def("set_shift", &T::set_shift)
        .def("set_id", &T::set_id)
        .def("set_mode", &T::set_mode)
        .def("set_scale", &T::set_scale)
        .def("set_bin", &T::set_bin)
        .def("set_owner", &T::set_owner);
}

namespace py = pybind11;

void bind_correlation_provider(py::module_& m) {
    auto corr = py::class_<CorrelationProvider, std::shared_ptr<CorrelationProvider>>(m, "CorrelationProvider");

    py::enum_<CorrelationProvider::CorrelationType>(corr, "CorrelationType")
        .value("STAT", CorrelationProvider::CorrelationType::STAT)
        .value("SYST", CorrelationProvider::CorrelationType::SYST)
        .value("COMBINED", CorrelationProvider::CorrelationType::COMBINED);

    corr
        .def(py::init<>())

        // ParamId / ParamId
        .def(
            "__call__",
            py::overload_cast<
                const ParamId&,
                const ParamId&,
                CorrelationProvider::CorrelationType
            >(&CorrelationProvider::operator()),
            py::arg("pid_1"),
            py::arg("pid_2"),
            py::arg("type")
        )

        // ExperimentObs / ExperimentObs
        .def(
            "__call__",
            py::overload_cast<
                const ExperimentObs&,
                const ExperimentObs&,
                CorrelationProvider::CorrelationType
            >(&CorrelationProvider::operator()),
            py::arg("obs_1"),
            py::arg("obs_2"),
            py::arg("type")
        )

        // experiment + Observables / Observables
        .def(
            "__call__",
            py::overload_cast<
                const std::string&,
                const Observables&,
                const Observables&,
                CorrelationProvider::CorrelationType
            >(&CorrelationProvider::operator()),
            py::arg("experiment"),
            py::arg("obs_1"),
            py::arg("obs_2"),
            py::arg("type")
        )

        // experiment + ObservableId / ObservableId
        .def(
            "__call__",
            py::overload_cast<
                const std::string&,
                const ObservableId&,
                const ObservableId&,
                CorrelationProvider::CorrelationType
            >(&CorrelationProvider::operator()),
            py::arg("experiment"),
            py::arg("obs_1"),
            py::arg("obs_2"),
            py::arg("type")
        )

        // experiment + BinnedObservableId / BinnedObservableId
        .def(
            "__call__",
            py::overload_cast<
                const std::string&,
                const BinnedObservableId&,
                const BinnedObservableId&,
                CorrelationProvider::CorrelationType
            >(&CorrelationProvider::operator()),
            py::arg("experiment"),
            py::arg("obs_1"),
            py::arg("obs_2"),
            py::arg("type")
        )

        // cross-experiment BinnedObservableId
        .def(
            "__call__",
            py::overload_cast<
                const std::string&,
                const BinnedObservableId&,
                const std::string&,
                const BinnedObservableId&,
                CorrelationProvider::CorrelationType
            >(&CorrelationProvider::operator()),
            py::arg("exp_1"),
            py::arg("obs_1"),
            py::arg("exp_2"),
            py::arg("obs_2"),
            py::arg("type")
        )

        // exists overloads
        .def(
            "exists",
            py::overload_cast<
                const ParamId&,
                const ParamId&,
                CorrelationProvider::CorrelationType
            >(&CorrelationProvider::exists),
            py::arg("pid_1"),
            py::arg("pid_2"),
            py::arg("type")
        )
        .def(
            "exists",
            py::overload_cast<
                const ExperimentObs&,
                const ExperimentObs&,
                CorrelationProvider::CorrelationType
            >(&CorrelationProvider::exists),
            py::arg("obs_1"),
            py::arg("obs_2"),
            py::arg("type")
        )
        .def(
            "exists",
            py::overload_cast<
                const std::string&,
                const Observables&,
                const Observables&,
                CorrelationProvider::CorrelationType
            >(&CorrelationProvider::exists),
            py::arg("experiment"),
            py::arg("obs_1"),
            py::arg("obs_2"),
            py::arg("type")
        )
        .def(
            "exists",
            py::overload_cast<
                const std::string&,
                const ObservableId&,
                const ObservableId&,
                CorrelationProvider::CorrelationType
            >(&CorrelationProvider::exists),
            py::arg("experiment"),
            py::arg("obs_1"),
            py::arg("obs_2"),
            py::arg("type")
        )
        .def(
            "exists",
            py::overload_cast<
                const std::string&,
                const BinnedObservableId&,
                const BinnedObservableId&,
                CorrelationProvider::CorrelationType
            >(&CorrelationProvider::exists),
            py::arg("experiment"),
            py::arg("obs_1"),
            py::arg("obs_2"),
            py::arg("type")
        )
        .def(
            "exists",
            py::overload_cast<
                const std::string&,
                const BinnedObservableId&,
                const std::string&,
                const BinnedObservableId&,
                CorrelationProvider::CorrelationType
            >(&CorrelationProvider::exists),
            py::arg("exp_1"),
            py::arg("obs_1"),
            py::arg("exp_2"),
            py::arg("obs_2"),
            py::arg("type")
        );
}

void init_core(py::module &m) {

    // py::class_<QCDHelper, std::shared_ptr<QCDHelper>>(m, "QCDHelper")
    //     .def_static(
    //         "mass_b_1S", &QCDHelper::mass_b_1S
    //     );

    //my technique to implicit cast from string to filesystem for python
    py::class_<std::filesystem::path>(m, "Path")
        .def(py::init<std::string>())
        .def("__str__", [](const std::filesystem::path &p) {
            return p.string();
        });

    py::implicitly_convertible<std::string, std::filesystem::path>();

    // ExternalFlag enum
    py::enum_<ExternalFlag>(m, "ExternalFlag")
    .value("IS_LHA_SPECTRUM", ExternalFlag::IS_LHA_SPECTRUM)
    .value("HAS_WILSON_INPUT", ExternalFlag::HAS_WILSON_INPUT)
    .value("HAS_TH_OBSERVABLE_INPUT", ExternalFlag::HAS_TH_OBSERVABLE_INPUT)
    .value("HYP_AS_SM_MARTY", ExternalFlag::HYP_AS_SM_MARTY)
    // .value("USE_MARTY", ExternalFlag::USE_MARTY)
    .export_values();

    // APIPath enum
    py::enum_<APIPath>(m, "APIPath")
    .value("LHA_PATH", APIPath::LHA_PATH)
    .export_values();
    
    py::enum_<ParameterMode>(m, "ParameterMode")
        .value("FIXED", ParameterMode::FIXED)
        .value("SHIFTABLE", ParameterMode::SHIFTABLE)
        .export_values();

    declare_parameter<Parameter>(m, "Parameter");

    py::class_<DependentParameter, Parameter, std::shared_ptr<DependentParameter>>(m, "DependentParameter")
        .def(py::init<ParamId, std::unordered_map<ParamId, std::shared_ptr<Parameter>>, DepParamUpdateFunc>())
        .def("init", &DependentParameter::init)
        .def("freeze", &DependentParameter::freeze)
        .def("unfreeze", &DependentParameter::unfreeze)
        .def("update", &DependentParameter::update)
        .def("dependsOn", &DependentParameter::dependsOn);

    // Config class
    py::class_<HyperisoConfig>(m, "HyperisoConfig")
    .def(py::init<>())
    .def_readwrite("flags", &HyperisoConfig::flags)
    .def_readwrite("model", &HyperisoConfig::model)
    .def_readwrite("mty_model_name", &HyperisoConfig::mty_model_name)
    .def_readwrite("mty_model_path", &HyperisoConfig::mty_model_path);

    // HyperisoMaster
    py::class_<HyperisoMaster, std::shared_ptr<HyperisoMaster>>(m, "HyperisoMaster")
    .def(py::init<>())
    .def("init", py::overload_cast<const std::string&, HyperisoConfig>(&HyperisoMaster::init))
    .def("init", py::overload_cast<const std::string&>(&HyperisoMaster::init))
    .def("check_flag", &HyperisoMaster::check_flag)
    .def("get_model", &HyperisoMaster::get_model)
    .def("switch_lha",  py::overload_cast<const std::string&, HyperisoConfig>(&HyperisoMaster::switch_lha));

    // ParameterSetter
    py::class_<ParameterShifter, std::shared_ptr<ParameterShifter>>(m, "ParameterShifter")
    .def(py::init<>())

    .def("mutate", &ParameterShifter::mutate, py::arg("pid"), py::arg("value"))
    .def("change_mode", &ParameterShifter::change_mode, py::arg("pid"), py::arg("mode"));

    py::class_<ParameterSetter, std::shared_ptr<ParameterSetter>>(m, "ParameterSetter")
    .def(py::init<>())

    .def("mutate", &ParameterSetter::mutate, py::arg("pid"), py::arg("value"))
    .def("change_mode", &ParameterSetter::change_mode, py::arg("pid"), py::arg("mode"));

    py::class_<ParameterProvider, std::shared_ptr<ParameterProvider>>(m, "ParameterProvider")
    .def(py::init<>())
    .def(py::init<ParameterType>(), py::arg("type"))

    .def("__call__",
         py::overload_cast<const ParamId&, DataType>(&ParameterProvider::operator(), py::const_),
         py::arg("pid"), py::arg("d_type") = DataType::VALUE)

    .def("__call__",
         py::overload_cast<const std::string&, const LhaID&, DataType>(&ParameterProvider::operator(), py::const_),
         py::arg("block"), py::arg("id"), py::arg("d_type") = DataType::VALUE)

    .def("exists",
         py::overload_cast<const ParamId&>(&ParameterProvider::exists, py::const_),
         py::arg("pid"))

    .def("exists",
         py::overload_cast<const std::string&, const LhaID&>(&ParameterProvider::exists, py::const_),
         py::arg("block"), py::arg("id"))

    .def("get_type", &ParameterProvider::get_type)

    .def("get_parameter", &ParameterProvider::get_parameter, py::arg("pid"))

    .def("__repr__", [](const ParameterProvider& self) {
        std::ostringstream oss;
        oss << "<ParameterProvider type=";
        if (auto t = self.get_type(); true) {
            oss << static_cast<int>(t);
        }
        oss << ">";
        return oss.str();
    });

    py::class_<BlockProvider>(m, "BlockProvider")
    .def(py::init<>())

    .def(
        "exists",
        &BlockProvider::exists,
        py::arg("blockname"),
        py::arg("type"),
        "Return whether a block exists for the given parameter type."
    )

    .def(
        "log_all_blocks",
        &BlockProvider::log_all_blocks,
        py::arg("type"),
        "Log all blocks for the given parameter type."
    )

    .def(
        "log_block",
        &BlockProvider::log_block,
        py::arg("type"),
        py::arg("blockname"),
        "Log one block for the given parameter type."
    )

    .def(
        "get_block",
        &BlockProvider::get_block,
        py::arg("type"),
        py::arg("blockname"),
        py::return_value_policy::move,
        "Return the content of one block for the given parameter type."
    )

    .def(
        "get_all_blocks",
        &BlockProvider::get_all_blocks,
        py::arg("type"),
        py::return_value_policy::move,
        "Return the names of all blocks for the given parameter type."
    );

    // py::enum_<CorrelationProvider::CorrelationType>(m, "CorrelationType")
    //     .value("STAT", CorrelationProvider::CorrelationType::STAT)
    //     .value("SYST", CorrelationProvider::CorrelationType::SYST)
    //     .value("COMBINED", CorrelationProvider::CorrelationType::COMBINED)
    //     .export_values();

    // // CorrelationProvider
    // py::class_<CorrelationProvider, std::shared_ptr<CorrelationProvider>>(m, "CorrelationProvider")
    // .def(py::init<>())
    // .def("correlation_from_paramid", py::overload_cast<const ParamId&, const ParamId&, CorrelationProvider::CorrelationType>(&CorrelationProvider::operator()), py::arg("pid_1"), py::arg("pid_2"), py::arg("type"))
    // .def("correlation_from_observable", py::overload_cast<const Observables&, const Observables&, CorrelationProvider::CorrelationType>(&CorrelationProvider::operator()), py::arg("pid_1"), py::arg("pid_2"), py::arg("type"));

    bind_correlation_provider(m);

    py::class_<DependantBlockInfoProvider, std::shared_ptr<DependantBlockInfoProvider>>(m, "DependantBlockInfoProvider")
        .def(py::init<>())
        .def(
            "is_dependent_block",
            &DependantBlockInfoProvider::is_dependent_block,
            py::arg("type"),
            py::arg("block_name"),
            "Return whether the requested block is implemented as a dependent block."
        )
        .def(
            "get_source_blocks",
            &DependantBlockInfoProvider::get_source_blocks,
            py::arg("type"),
            py::arg("block_name"),
            "Return the direct upstream/source blocks used by the requested block."
        )
        .def(
            "get_dependent_blocks",
            &DependantBlockInfoProvider::get_dependent_blocks,
            py::arg("type"),
            py::arg("block_name"),
            "Return the direct downstream blocks depending on the requested block."
        )
        .def(
            "get_all_source_blocks",
            &DependantBlockInfoProvider::get_all_source_blocks,
            py::arg("type"),
            py::arg("block_name"),
            "Return all transitive upstream/source blocks used by the requested block."
        )
        .def(
            "get_all_dependent_blocks",
            &DependantBlockInfoProvider::get_all_dependent_blocks,
            py::arg("type"),
            py::arg("block_name"),
            "Return all transitive downstream blocks depending on the requested block."
        );

    // QCDProvider
    py::class_<QCDConstants>(m, "QCDConstants")
        .def_readonly_static("Nc", &QCDConstants::Nc)
        .def_readonly_static("C_F", &QCDConstants::C_F)
        .def_readonly_static("C_A", &QCDConstants::C_A)
        .def_readonly_static("beta", &QCDConstants::beta)
        .def_readonly_static("gamma", &QCDConstants::gamma);

    
        
    py::class_<QCDProvider, std::shared_ptr<QCDProvider>>(m, "QCDProvider")
        .def(py::init<>())
        .def("compute_alphas", py::overload_cast<AlphasConfig>(&QCDProvider::operator()), py::arg("alpha_config"))
        .def("compute_mass", py::overload_cast<MassConfig>(&QCDProvider::operator()), py::arg("mass_config"))
        .def("get_constants", [](QCDProvider &self) {
            return self.get_constants();
        }, py::return_value_policy::reference);
    
    py::class_<QEDProvider, std::shared_ptr<QEDProvider>>(m, "QEDProvider")
        .def(py::init<>())
        .def("compute_alphaem", py::overload_cast<AlphasConfig>(&QEDProvider::operator()), py::arg("alpha_config"));
        
    // APIAdapter
    py::class_<APIAdapter, std::shared_ptr<APIAdapter>>(m, "APIAdapter")
    .def(py::init<>())
    .def("check_flag", &APIAdapter::check_flag)
    .def("get_path", &APIAdapter::get_path)
    .def("get_all_blocks", &APIAdapter::get_all_blocks)
    .def("get_blocks_list", &APIAdapter::get_blocks_list, py::arg("param_type") = ParameterType::SM)
    .def("get_block_infos", &APIAdapter::get_block_infos, py::arg("block"), py::arg("param_type") = ParameterType::SM)
    .def("get_type_of_block", &APIAdapter::get_type_of_block);


    py::class_<Freezer>(m, "Freezer")
    .def_static("freeze_block", 
                static_cast<void(*)(const ParameterType&, const std::string&)>(&Freezer::freeze),
                py::arg("p_type"), py::arg("block_name"))
    .def_static("freeze_param", 
                static_cast<void(*)(const ParamId&)>(&Freezer::freeze),
                py::arg("pid"))
    .def_static("unfreeze_block", 
                static_cast<void(*)(const ParameterType&, const std::string&)>(&Freezer::unfreeze),
                py::arg("p_type"), py::arg("block_name"))
    .def_static("unfreeze_param", 
                static_cast<void(*)(const ParamId&)>(&Freezer::unfreeze),
                py::arg("pid"));

    // DBMemento
    py::class_<DBMemento, std::shared_ptr<DBMemento>>(m, "DBMemento")
        .def(py::init<>())
        .def("takeSnapshot", &DBMemento::takeSnapshot, py::arg("blocks"))
        .def("restore", &DBMemento::restore, py::arg("n_steps") = 1)
        .def("print_snapshot_content", &DBMemento::print_snapshot_content, py::arg("n") = 1)
        .def("stack_size", &DBMemento::stack_size);

            py::enum_<MartyPath>(m, "MartyPath")
        .value("MODEL_FILE", MartyPath::MODEL_FILE)
        .value("TEMPLATE_DIR", MartyPath::TEMPLATE_DIR)
        .value("PARAM_MAPPING_DIR", MartyPath::PARAM_MAPPING_DIR)
        .export_values();
    

    py::class_<MartyAdapter, std::shared_ptr<MartyAdapter>>(m, "MartyAdapter")
        .def(py::init<>())

        .def("get_path", [](MartyAdapter& self, MartyPath path) {
            return self.get_path(path).string();
        }, py::arg("path"))

        .def("check_flag", &MartyAdapter::check_flag, py::arg("flag"))

        .def("get_marty_model_name", &MartyAdapter::get_marty_model_name)

        .def("__repr__", [](const MartyAdapter& self) {
            return "<MartyAdapter: model=" + self.get_marty_model_name() + ">";
        });

}

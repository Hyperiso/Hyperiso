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
#include "DependencyPruner.h"
#include "FileWriter.h"

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
    .value("ASSETS_ROOT", APIPath::ASSETS_ROOT)
    .value("DEFAULT_PARAM_VALUES", APIPath::DEFAULT_PARAM_VALUES)
    .value("DEFAULT_OBS_VALUES", APIPath::DEFAULT_OBS_VALUES)
    .value("DEFAULT_PARAM_CORR", APIPath::DEFAULT_PARAM_CORR)
    .value("DEFAULT_OBS_CORR", APIPath::DEFAULT_OBS_CORR)
    .value("DEFAULT_NUISANCES", APIPath::DEFAULT_NUISANCES)
    .value("USER_SM_PARAMS", APIPath::USER_SM_PARAMS)
    .value("USER_FLAVOR_PARAMS", APIPath::USER_FLAVOR_PARAMS)
    .value("USER_DECAY_PARAMS", APIPath::USER_DECAY_PARAMS)
    .value("USER_OBS_VALUES", APIPath::USER_OBS_VALUES)
    .value("USER_PARAM_CORR", APIPath::USER_PARAM_CORR)
    .value("USER_OBS_CORR", APIPath::USER_OBS_CORR)
    .value("USER_NUISANCES", APIPath::USER_NUISANCES)
    .value("PARAM_MAPPING_DIR", APIPath::PARAM_MAPPING_DIR)
    .value("TEMPLATE_DIR", APIPath::TEMPLATE_DIR)
    .value("SPECTRUM_DIR", APIPath::SPECTRUM_DIR)
    .value("MARTY_TEMP_DIR", APIPath::MARTY_TEMP_DIR)
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
    .def_readwrite("mty_model_path", &HyperisoConfig::mty_model_path)
    .def_readwrite("mty_bsm_mapping_path", &HyperisoConfig::mty_bsm_mapping_path)
    .def_readwrite("mty_order_policy", &HyperisoConfig::mty_order_policy)
    .def_readwrite("mty_tree_fermion_order", &HyperisoConfig::mty_tree_fermion_order);

    // HyperisoMaster
    py::class_<HyperisoMaster, std::shared_ptr<HyperisoMaster>>(m, "HyperisoMaster")
        .def(py::init<>())
        .def("init", py::overload_cast<const std::string&, HyperisoConfig>(&HyperisoMaster::init))
        .def("init", py::overload_cast<const std::string&>(&HyperisoMaster::init))
        .def(
            "pre_init_add_block",
            [](HyperisoMaster& self,
            const std::string& block_name,
            size_t item_count,
            size_t value_idx,
            int scale_idx,
            int rg_idx,
            int bin_idx,
            bool global_scale) {
                self.pre_init_add_block(BlockName(block_name),
                                        item_count,
                                        value_idx,
                                        scale_idx,
                                        rg_idx,
                                        bin_idx,
                                        global_scale);
            },
            py::arg("block_name"),
            py::arg("item_count") = 2,
            py::arg("value_idx") = 1,
            py::arg("scale_idx") = -1,
            py::arg("rg_idx") = -1,
            py::arg("bin_idx") = -1,
            py::arg("global_scale") = false,
            R"doc(Register an additional LHA block prototype before initialization.

    The block is added to the LHA parser prototype registry and is routed to the
    BSM parameter namespace by default. Call this before init(). Multiple custom
    blocks can be registered by calling this method several times.)doc"
        )
        .def(
            "pre_init_set_marty_path",
            &HyperisoMaster::pre_init_set_marty_path,
            py::arg("marty_path"),
            R"doc(Register an existing MARTY installation before initialization.

    The path may point to the MARTY install prefix, to a directory containing
    MARTY_INSTALL/ or install/, to include/, to lib/, to marty.h, or to a libmarty
    library file. Hyperiso validates include/marty.h and lib/libmarty.* before MARTY
    mode is allowed to initialize.)doc"
        )
        .def(
            "pre_init_set_softsusy_path",
            &HyperisoMaster::pre_init_set_softsusy_path,
            py::arg("softsusy_path"),
            R"doc(Register a SOFTSUSY executable or installation directory before initialization.

    The path may point directly to softpoint.x, or to a directory containing
    softpoint.x, bin/softpoint.x, or src/SOFTSUSY/softpoint.x. This is the
    recommended Python path when Hyperiso was not built with bundled
    SOFTSUSY support.)doc"
        )
        .def(
            "pre_init_set_paths",
            [](HyperisoMaster& self, const std::map<APIPath, std::string>& path_overrides) {
                self.pre_init_set_paths(path_overrides);
            },
            py::arg("path_overrides"),
            R"doc(Override selected Hyperiso filesystem paths before initialization.

    Args:
        path_overrides: Mapping from APIPath to filesystem path. Default input
            files must be existing .json files, user input files must be existing
            .yaml/.yml files, and directory entries must be existing directories.

    LHA_PATH is intentionally not accepted here because the active LHA file is
    provided through init() or switch_lha().)doc"
        )
        .def(
            "pre_init_set_marty_cache_dir",
            &HyperisoMaster::pre_init_set_marty_cache_dir,
            py::arg("cache_dir"),
            R"doc(Set the writable MARTY generated-code/cache directory before initialization.

    The directory is created by the C++ layer if it does not exist.)doc"
        )
        .def(
            "pre_init_set_spectrum_cache_dir",
            &HyperisoMaster::pre_init_set_spectrum_cache_dir,
            py::arg("cache_dir"),
            R"doc(Set the writable spectrum cache directory before initialization.

    The directory is created by the C++ layer if it does not exist.)doc"
        )
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

    py::class_<FileWriter>(
        m,
        "DatabaseWriter",
        R"doc(Export the initialized Core parameter database.

The output format is selected from the destination suffix. Supported suffixes
are .json, .yaml, .yml, .lha, .slha and .flha. Exports can contain the whole
database, selected blocks, or selected block/id parameters.)doc"
    )
        .def(py::init<>())
        .def(
            "write",
            [](FileWriter& writer, const std::string& destination) {
                writer.write(destination);
            },
            py::arg("destination"),
            R"doc(Export the complete initialized database.

Args:
    destination: Output path. The filename suffix selects JSON, YAML or LHA.)doc"
        )
        .def(
            "write_blocks",
            [](FileWriter& writer,
               const std::string& destination,
               const std::vector<std::string>& block_names) {
                std::vector<BlockName> names;
                names.reserve(block_names.size());
                for (const auto& name : block_names) {
                    names.emplace_back(name);
                }
                writer.write_blocks(destination, names);
            },
            py::arg("destination"),
            py::arg("block_names"),
            R"doc(Export selected blocks from the initialized database.

Args:
    destination: Output path. The filename suffix selects the format.
    block_names: Block names or aliases to export.)doc"
        )
        .def(
            "write_parameters",
            [](FileWriter& writer,
               const std::string& destination,
               const std::vector<ParamId>& parameter_ids) {
                writer.write_parameters(destination, parameter_ids);
            },
            py::arg("destination"),
            py::arg("parameter_ids"),
            R"doc(Export selected parameters addressed by block and LHA id.

Args:
    destination: Output path. The filename suffix selects the format.
    parameter_ids: ParamId objects identifying the entries to export.)doc"
        );

    // Backward/adapter-oriented name for callers working directly with the
    // C++ abstraction. The Python wrapper uses the clearer DatabaseWriter name.
    m.attr("FileWriter") = m.attr("DatabaseWriter");

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
    
    py::class_<DependencyPruner, std::shared_ptr<DependencyPruner>>(
        m,
        "DependencyPruner",
        R"doc(Adapter used to detach and reattach dependency links in the parameter graph.

This binding exposes the C++ DependencyPruner adapter to Python. Block names are
accepted as strings and converted to BlockName on the C++ side. Parameter ids are
expected to be bound LhaID instances.)doc"
    )
        .def(
            py::init<>(),
            R"doc(Create a dependency-pruning adapter backed by Parameters.)doc"
        )
        .def(
            "detach_block",
            [](DependencyPruner& self, ParameterType type, const std::string& block_name) {
                self.detach_block(type, BlockName(block_name));
            },
            py::arg("type"),
            py::arg("block_name"),
            R"doc(Detach a dependent block from its upstream source blocks.

Args:
    type: Parameter namespace containing the block.
    block_name: Name of the dependent block to detach.)doc"
        )
        .def(
            "reattach_block",
            [](DependencyPruner& self, ParameterType type, const std::string& block_name) {
                self.reattach_block(type, BlockName(block_name));
            },
            py::arg("type"),
            py::arg("block_name"),
            R"doc(Reattach a previously detached dependent block to its sources.

Args:
    type: Parameter namespace containing the block.
    block_name: Name of the dependent block to reattach.)doc"
        )
        .def(
            "detach_parameter",
            [](DependencyPruner& self, ParameterType type, const std::string& block_name, const LhaID& id) {
                self.detach_parameter(type, BlockName(block_name), id);
            },
            py::arg("type"),
            py::arg("block_name"),
            py::arg("id"),
            R"doc(Detach a dependent parameter from its upstream source parameters.

Args:
    type: Parameter namespace containing the parameter.
    block_name: Name of the block containing the parameter.
    id: LHA identifier of the dependent parameter.)doc"
        )
        .def(
            "reattach_parameter",
            [](DependencyPruner& self, ParameterType type, const std::string& block_name, const LhaID& id) {
                self.reattach_parameter(type, BlockName(block_name), id);
            },
            py::arg("type"),
            py::arg("block_name"),
            py::arg("id"),
            R"doc(Reattach a previously detached dependent parameter to its sources.

Args:
    type: Parameter namespace containing the parameter.
    block_name: Name of the block containing the parameter.
    id: LHA identifier of the dependent parameter.)doc"
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
        .value("SM_MODEL_FILE", MartyPath::SM_MODEL_FILE)
        .value("TEMPLATE_DIR", MartyPath::TEMPLATE_DIR)
        .value("PARAM_MAPPING_DIR", MartyPath::PARAM_MAPPING_DIR)
        .value("SM_MAPPING_FILE", MartyPath::SM_MAPPING_FILE)
        .value("BSM_MAPPING_FILE", MartyPath::BSM_MAPPING_FILE)
        .value("MARTY_TEMP_DIR", MartyPath::MARTY_TEMP_DIR)
        .export_values();
    

    py::class_<MartyAdapter, std::shared_ptr<MartyAdapter>>(m, "MartyAdapter")
        .def(py::init<>())

        .def("get_path", [](MartyAdapter& self, MartyPath path) {
            return self.get_path(path).string();
        }, py::arg("path"))

        .def("get_optional_path", [](MartyAdapter& self, MartyPath path) -> py::object {
            auto resolved = self.get_optional_path(path);
            if (!resolved.has_value()) {
                return py::none();
            }
            return py::str(resolved->string());
        }, py::arg("path"))

        .def("check_flag", &MartyAdapter::check_flag, py::arg("flag"))

        .def("get_marty_model_name", &MartyAdapter::get_marty_model_name)
        .def("get_marty_order_policy", &MartyAdapter::get_marty_order_policy)

        .def("__repr__", [](const MartyAdapter& self) {
            return "<MartyAdapter: model=" + self.get_marty_model_name() + ">";
        });

}

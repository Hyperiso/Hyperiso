#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <optional>
#include <filesystem>
#include "General.h"
// #include "MemoryManager.h"
// #include "Parameters.h"
#include "HyperisoMaster.h"
#include "Config.h"
#include "ParameterSetter.h"
#include "ParameterProvider.h"
#include "CorrelationProvider.h"
#include "QCDProvider.h"
#include "Freezer.h"
#include "DBMemento.h"
#include "APIAdapter.h"

namespace py = pybind11;

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
    .value("USE_MARTY", ExternalFlag::USE_MARTY)
    .export_values();

    // APIPath enum
    py::enum_<APIPath>(m, "APIPath")
    .value("LHA_PATH", APIPath::LHA_PATH)
    .export_values();

    // Config class
    py::class_<Config>(m, "Config")
    .def(py::init<>())
    .def_readwrite("flags", &Config::flags)
    .def_readwrite("model", &Config::model)
    .def_readwrite("mty_model_name", &Config::mty_model_name)
    .def_readwrite("mty_model_path", &Config::mty_model_path);

    // HyperisoMaster
    py::class_<HyperisoMaster, std::shared_ptr<HyperisoMaster>>(m, "HyperisoMaster")
    .def(py::init<>())
    .def("init", py::overload_cast<const std::string&, Config>(&HyperisoMaster::init))
    .def("init", py::overload_cast<const std::string&>(&HyperisoMaster::init))
    .def("check_flag", &HyperisoMaster::check_flag)
    .def("get_model", &HyperisoMaster::get_model);

    // ParameterSetter
    py::class_<ParameterSetter, std::shared_ptr<ParameterSetter>>(m, "ParameterSetter")
    .def(py::init<>())
    .def("mutate", &ParameterSetter::mutate);

    // ParameterProvider::DataType
    py::enum_<ParameterProvider::DataType>(m, "DataType")
    .value("VALUE", ParameterProvider::DataType::VALUE)
    .value("STD_STAT", ParameterProvider::DataType::STD_STAT)
    .value("STD_SYST", ParameterProvider::DataType::STD_SYST)
    .value("STD_COMBINED", ParameterProvider::DataType::STD_COMBINED)
    .export_values();

    // ParameterProvider
    py::class_<ParameterProvider, std::shared_ptr<ParameterProvider>>(m, "ParameterProvider")
    .def(py::init<>())
    .def(py::init<ParameterType>())
    .def("__call__", py::overload_cast<const ParamId&, ParameterProvider::DataType>(&ParameterProvider::operator()), py::arg("pid"), py::arg("d_type") = ParameterProvider::DataType::VALUE)
    .def("__call__", py::overload_cast<const std::string&, const LhaID&, ParameterProvider::DataType>(&ParameterProvider::operator(), py::const_), py::arg("block"), py::arg("id"), py::arg("d_type") = ParameterProvider::DataType::VALUE)
    .def("exists", py::overload_cast<const ParamId&>(&ParameterProvider::exists, py::const_))
    .def("exists", py::overload_cast<const std::string&, const LhaID&>(&ParameterProvider::exists, py::const_))
    .def("get_type", &ParameterProvider::get_type);


    // CorrelationProvider
    py::class_<CorrelationProvider, std::shared_ptr<CorrelationProvider>>(m, "CorrelationProvider")
    .def(py::init<>())
    .def("__call__", py::overload_cast<const ParamId&, const ParamId&, CorrelationProvider::CorrelationType>(&CorrelationProvider::operator()), py::arg("pid_1"), py::arg("pid_2"), py::arg("type"))
    .def("__call__", py::overload_cast<const Observables&, const Observables&, CorrelationProvider::CorrelationType>(&CorrelationProvider::operator()), py::arg("pid_1"), py::arg("pid_2"), py::arg("type"));

    // QCDProvider
    py::class_<QCDProvider, std::shared_ptr<QCDProvider>>(m, "QCDProvider")
    .def(py::init<>())
    .def("__call__", py::overload_cast<AlphasConfig>(&QCDProvider::operator()), py::arg("alpha_config"))
    .def("__call__", py::overload_cast<MassConfig>(&QCDProvider::operator()), py::arg("mass_config"))
    .def("get_constants", &QCDProvider::get_constants);

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

}

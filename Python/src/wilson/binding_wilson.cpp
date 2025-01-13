#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#include "WilsonManager.h"
#include "MartyWilson.h"
#include "WilsonInterface.h"
#include "Wilson_susyv2.h"
#include "Wilson_THDMv2.h"
namespace py = pybind11;

// Initialisation des groupes de coefficients
void init_coefficient_groups(py::module &m) {
    py::class_<CoefficientGroup, std::shared_ptr<CoefficientGroup>>(m, "CoefficientGroup")
        .def("get_matching", &CoefficientGroup::getMatching)
        .def("get_run", &CoefficientGroup::getRun)
        .def("set_q_match", &CoefficientGroup::set_Q_match)
        .def("set_q_run", &CoefficientGroup::set_Q_run);

    py::class_<BCoefficientGroup, CoefficientGroup, std::shared_ptr<BCoefficientGroup>>(m, "BCoefficientGroup")
        .def(py::init<>());

    py::class_<BCoefficientGroupMarty, CoefficientGroup, std::shared_ptr<BCoefficientGroupMarty>>(m, "BCoefficientGroupMarty")
        .def(py::init<>());

    py::class_<BPrimeCoefficientGroup, CoefficientGroup, std::shared_ptr<BPrimeCoefficientGroup>>(m, "BPrimeCoefficientGroup")
        .def(py::init<>());

    py::class_<BPrimeCoefficientGroupMarty, CoefficientGroup, std::shared_ptr<BPrimeCoefficientGroupMarty>>(m, "BPrimeCoefficientGroupMarty")
        .def(py::init<>());

    py::class_<BScalarCoefficientGroup, CoefficientGroup, std::shared_ptr<BScalarCoefficientGroup>>(m, "BScalarCoefficientGroup")
        .def(py::init<>());

    py::class_<BScalarCoefficientGroupMarty, CoefficientGroup, std::shared_ptr<BScalarCoefficientGroupMarty>>(m, "BScalarCoefficientGroupMarty")
        .def(py::init<>());

    py::class_<BCoefficientGroup_susy, CoefficientGroup, std::shared_ptr<BCoefficientGroup_susy>>(m, "BCoefficientGroup_susy")
        .def(py::init<>());

    py::class_<BCoefficientGroup_THDM, CoefficientGroup, std::shared_ptr<BCoefficientGroup_THDM>>(m, "BCoefficientGroup_THDM")
        .def(py::init<>());

    py::class_<BPrimeCoefficientGroup_susy, CoefficientGroup, std::shared_ptr<BPrimeCoefficientGroup_susy>>(m, "BPrimeCoefficientGroup_susy")
        .def(py::init<>());

    py::class_<BPrimeCoefficientGroup_THDM, CoefficientGroup, std::shared_ptr<BPrimeCoefficientGroup_THDM>>(m, "BPrimeCoefficientGroup_THDM")
        .def(py::init<>());

    py::class_<BScalarCoefficientGroup_susy, CoefficientGroup, std::shared_ptr<BScalarCoefficientGroup_susy>>(m, "BScalarCoefficientGroup_susy")
        .def(py::init<>());

    py::class_<BScalarCoefficientGroup_THDM, CoefficientGroup, std::shared_ptr<BScalarCoefficientGroup_THDM>>(m, "BScalarCoefficientGroup_THDM")
        .def(py::init<>());
}

// Initialisation des coefficients Wilson
void init_wilson_coefficient(py::module &m) {
    py::class_<WilsonCoefficient, std::shared_ptr<WilsonCoefficient>>(m, "WilsonCoefficient")
        .def("get_coefficient_matching_value", &WilsonCoefficient::get_CoefficientMatchingValue)
        .def("get_coefficient_run_value", &WilsonCoefficient::get_CoefficientRunValue)
        .def("set_coefficient_matching_value", &WilsonCoefficient::set_CoefficientMatchingValue)
        .def("set_coefficient_run_value", &WilsonCoefficient::set_WilsonCoeffRun)
        .def("get_q_match", &WilsonCoefficient::get_Q_match)
        .def("get_q", &WilsonCoefficient::get_Q);
}

// Initialisation du gestionnaire de coefficients
void init_coefficient_manager(py::module &m) {
    // py::enum_<Model>(m, "Model")
    // .value("SM", Model::SM)
    // .value("THDM", Model::THDM)
    // .value("SUSY", Model::THDM)
    // .value("CUSTOM", Model::THDM)
    // .export_values();

    py::class_<CoefficientManager, std::shared_ptr<CoefficientManager>>(m, "CoefficientManager")
        .def_static("get_instance", [](const std::string &modelName) {
            return CoefficientManager::GetInstance(modelName);
        }, py::return_value_policy::reference)

        .def_static("initialize", &CoefficientManager::initialize, py::arg("lhaFile"), py::arg_v("model",Model::SM),
            py::arg_v("use_marty",false),
            py::arg_v("is_spectrum", false),
            py::arg_v("has_wilsons", false),
            py::arg_v("has_obs", false),
            "Initialise MemoryManager avec un fichier LHA et un modèle.")
        .def("register_coefficient_group", &CoefficientManager::registerCoefficientGroup)
        .def("get_state", &CoefficientManager::get_state)
        .def("get_alpha_s", &CoefficientManager::getAlphaS)
        .def("set_q_match", &CoefficientManager::setQMatch)
        .def("set_params", &CoefficientManager::setParams)
        .def("get_params", &CoefficientManager::get_params)
        .def("set_group_scale", &CoefficientManager::setGroupScale)
        .def("set_matching_coefficient", &CoefficientManager::setMatchingCoefficient)
        .def("set_run_coefficient", &CoefficientManager::setRunCoefficient)
        .def("get_matching_coefficient", &CoefficientManager::getMatchingCoefficient)
        .def("get_run_coefficient", &CoefficientManager::getRunCoefficient)
        .def("get_coefficient_group", &CoefficientManager::getCoefficientGroup, py::return_value_policy::reference);
}

// Initialisation des paramètres Wilson
void init_wilson_parameters(py::module &m) {
    py::class_<Wilson_parameters, std::shared_ptr<Wilson_parameters>>(m, "WilsonParameters")
        .def_static("get_instance", &Wilson_parameters::GetInstance, py::return_value_policy::reference)
        .def("set_mu", &Wilson_parameters::SetMu)
        .def("set_mu_w", &Wilson_parameters::SetMuW)
        .def("set_gen", &Wilson_parameters::set_gen);
}

void init_wilson_interface(py::module &m) {
    py::enum_<BWilsonCoefficients>(m, "BWilsonCoefficients")
        .value("C1", BWilsonCoefficients::C1)
        .value("C2", BWilsonCoefficients::C2)
        .value("C3", BWilsonCoefficients::C3)
        .value("C4", BWilsonCoefficients::C4)
        .value("C5", BWilsonCoefficients::C5)
        .value("C6", BWilsonCoefficients::C6)
        .value("C7", BWilsonCoefficients::C7)
        .value("C8", BWilsonCoefficients::C8)
        .value("C9", BWilsonCoefficients::C9)
        .value("C10", BWilsonCoefficients::C10)
        .value("CQ1", BWilsonCoefficients::CQ1)
        .value("CQ2", BWilsonCoefficients::CQ2)
        .value("CP1", BWilsonCoefficients::CP1)
        .value("CP2", BWilsonCoefficients::CP2)
        .value("CP3", BWilsonCoefficients::CP3)
        .value("CP4", BWilsonCoefficients::CP4)
        .value("CP5", BWilsonCoefficients::CP5)
        .value("CP6", BWilsonCoefficients::CP6)
        .value("CP7", BWilsonCoefficients::CP7)
        .value("CP8", BWilsonCoefficients::CP8)
        .value("CP9", BWilsonCoefficients::C9)
        .value("CP10", BWilsonCoefficients::C10)
        .value("CPQ1", BWilsonCoefficients::CPQ1)
        .value("CPQ2", BWilsonCoefficients::CPQ2)
        .export_values();

    py::enum_<WilsonGroups>(m, "WilsonGroups")
        .value("BCoefficients", WilsonGroups::BCoefficients)
        .value("BPrimeCoefficients", WilsonGroups::BPrimeCoefficients)
        .value("BScalarCoefficients", WilsonGroups::BScalarCoefficients)
        .value("BCoefficients_THDM", WilsonGroups::BCoefficients_THDM)
        .value("BPrimeCoefficients_THDM", WilsonGroups::BPrimeCoefficients_THDM)
        .value("BScalarCoefficients_THDM", WilsonGroups::BScalarCoefficients_THDM)
        .value("BCoefficients_SUSY", WilsonGroups::BCoefficients_SUSY)
        .value("BPrimeCoefficients_SUSY", WilsonGroups::BPrimeCoefficients_SUSY)
        .value("BScalarCoefficients_SUSY", WilsonGroups::BScalarCoefficients_SUSY)
        .value("BCoefficients_MARTY", WilsonGroups::BCoefficients_MARTY)
        .value("BPrimeCoefficients_MARTY", WilsonGroups::BPrimeCoefficients_MARTY)
        .value("BScalarCoefficients_MARTY", WilsonGroups::BScalarCoefficients_MARTY)
        .export_values();

    py::enum_<BWilsonBasis>(m, "BWilsonBasis")
        .value("STANDARD", BWilsonBasis::STANDARD)
        .value("TRADITIONAL", BWilsonBasis::TRADITIONAL)
        .export_values();

    py::class_<WilsonInterface, std::shared_ptr<WilsonInterface>>(m, "WilsonInterface")
        .def(py::init<const std::string &>())
        .def("add_wilson_group", &WilsonInterface::AddWilsonGroup)
        .def("set_q_match", &WilsonInterface::setQMatch)
        .def("set_params", &WilsonInterface::setParams)
        .def("set_matching_coefficient", &WilsonInterface::setMatchingCoefficient)
        .def("set_group_scale", &WilsonInterface::setGroupScale)
        .def("set_run_coefficient", &WilsonInterface::setRunCoefficient)
        .def("switch_basis", &WilsonInterface::switchbasis)
        .def("get_alpha_s", &WilsonInterface::getAlphaS)
        .def("get_matching_coefficient", &WilsonInterface::getMatchingCoefficient)
        .def("builder", &WilsonInterface::Builder);
}

// Fonction principale d'initialisation pour le module Wilson
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

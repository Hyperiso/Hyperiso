#include "MartyFileWriter.h"
#include "../MartyNumericalPolicy.h"

#include <iomanip>
#include <sstream>
#include <unordered_set>

namespace {
std::string cpp_double_literal(double value) {
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

const std::unordered_set<std::string>& wilsons_without_mudim() {
    static const std::unordered_set<std::string> values = {
        "C10", "CP10"
    };
    return values;
}

const std::unordered_set<std::string>& wilsons_with_marty_split_sm_components() {
    static const std::unordered_set<std::string> values = {
        // Empty by design: SM primed semileptonic coefficients are supplied by
        // the builtin backend.  For C9/CP9/CP10 the MARTY split is BSM-only.
    };
    return values;
}

bool should_read_marty_split_sm_components(const std::string& wilson) {
    return wilsons_with_marty_split_sm_components().find(wilson) != wilsons_with_marty_split_sm_components().end();
}

bool uses_photon_veto_diagnostic_policy(const std::string& wilson) {
    return wilson == "C9" || wilson == "CP9";
}

bool uses_split_regprop_policy(const std::string& wilson) {
    // CP10 still uses the newer two-piece split.  C9/CP9 deliberately keep the
    // older stable policy: the raw photon-linker projection is diagnostic only
    // and is never added to the physical Wilson coefficient.
    return wilson == "CP10";
}
}

MartyFileWriter::MartyFileWriter(const std::string& wilson,
                                 const std::string& model,
                                 bool bsm_split_generation,
                                 bool full_target_generation) :
    wilson(wilson),
    model(model),
    bsm_split_generation(bsm_split_generation),
    full_target_generation(full_target_generation) {}

bool MartyFileWriter::should_set_mudim() const {
    return wilsons_without_mudim().find(this->wilson) == wilsons_without_mudim().end();
}

void MartyFileWriter::add_output_writer(std::ofstream& outputFile) {
    const std::string default_regprop = cpp_double_literal(MartyNumericalPolicy::kDefaultRegProp);
    const std::string photon_regprop = cpp_double_literal(MartyNumericalPolicy::kPhotonDiagnosticRegProp);

    outputFile << "\tconst std::string& path = output_file_path;\n";

    if (should_set_mudim()) {
        outputFile << "\tsetMu(Q_match);\n";
    } else {
        outputFile << "\t// setMu(Q_match) intentionally skipped for " << this->wilson
                   << " because LoopTools mudim must stay at its default value.\n";
    }

    if (bsm_split_generation && uses_photon_veto_diagnostic_policy(this->wilson)) {
        outputFile << "\t// Stable C9/CP9 MARTY policy: ordinary/non-photon BSM matching uses the small regulator.\n";
        outputFile << "\t// The raw massless-photon linker is disabled by default; --raw-photon-diagnostic uses the dedicated photon regulator.\n";
        outputFile << "\t// A finite model-specific photon term, when known, is added later through WilsonMatchingPatch.\n";
        outputFile << "\tauto hyperiso_regprop_it = param.realParams.find(\"reg_prop\");\n";
        outputFile << "\tbool hyperiso_has_regprop = hyperiso_regprop_it != param.realParams.end() && hyperiso_regprop_it->second != nullptr;\n";
        outputFile << "\tdouble hyperiso_regprop_saved = " << default_regprop << ";\n";
        outputFile << "\tif (hyperiso_has_regprop) {\n";
        outputFile << "\t\thyperiso_regprop_saved = static_cast<double>(*hyperiso_regprop_it->second);\n";
        outputFile << "\t\t*hyperiso_regprop_it->second = " << default_regprop << ";\n";
        outputFile << "\t}\n";
        outputFile << "\tauto hyperiso_bsm_non_photon = " + wilson + "(param);\n";
        outputFile << "\tauto hyperiso_bsm_photon_raw = 0.0 * hyperiso_bsm_non_photon;\n";
        outputFile << "\tif (raw_photon_diagnostic) {\n";
        outputFile << "\t\tif (hyperiso_has_regprop) {\n";
        outputFile << "\t\t\t*hyperiso_regprop_it->second = " << photon_regprop << ";\n";
        outputFile << "\t\t}\n";
        outputFile << "\t\thyperiso_bsm_photon_raw = " + wilson + "_A(param);\n";
        outputFile << "\t}\n";
        outputFile << "\tif (hyperiso_has_regprop) {\n";
        outputFile << "\t\t*hyperiso_regprop_it->second = hyperiso_regprop_saved;\n";
        outputFile << "\t}\n";
        outputFile << "\tauto hyperiso_bsm_physical = hyperiso_bsm_non_photon;\n";
        outputFile << "\tauto hyperiso_raw_photon_sum = hyperiso_bsm_non_photon + hyperiso_bsm_photon_raw;\n";
        outputFile << "\tauto hyperiso_zero = 0.0 * hyperiso_bsm_non_photon;\n";
        if (full_target_generation) {
            outputFile << "\twriteWilsonCoefficients(\"" + wilson + "\", hyperiso_bsm_physical, Q_match, path);\n";
            outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_TARGET_SPLIT\", hyperiso_bsm_physical, Q_match, path);\n";
        } else {
            outputFile << "\twriteWilsonCoefficients(\"" + wilson + "\", hyperiso_bsm_physical, Q_match, path);\n";
        }
        // Keep the historical diagnostic names, but make all names advertised as
        // physical/split point to the photon-vetoed value.  The unsafe raw sum is
        // exported under an explicit RAW name so it cannot be consumed by mistake.
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_BSM_SPLIT\", hyperiso_bsm_physical, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_SM_SPLIT\", hyperiso_zero, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_TOTAL_SPLIT\", hyperiso_bsm_physical, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_NONPHOTON\", hyperiso_bsm_non_photon, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_A\", hyperiso_bsm_photon_raw, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_RAW_PHOTON_SUM\", hyperiso_raw_photon_sum, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_SM_NONPHOTON\", hyperiso_zero, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_SM_A\", hyperiso_zero, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_VECTOR\", hyperiso_bsm_non_photon, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_SCALAR\", hyperiso_zero, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_SM_VECTOR\", hyperiso_zero, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_SM_SCALAR\", hyperiso_zero, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_TOTAL_VECTOR\", hyperiso_bsm_non_photon, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_TOTAL_SCALAR\", hyperiso_zero, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_SM_COMPONENT\", hyperiso_zero, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_TOTAL_COMPONENT\", hyperiso_bsm_physical, Q_match, path);\n";
        return;
    }

    if (bsm_split_generation && uses_split_regprop_policy(this->wilson)) {
        outputFile << "\t// Split MARTY coefficient policy. " << wilson
                   << " is generated as non-photon and photon-linker pieces.\n";
        outputFile << "\t// Non-photon and photon-linker pieces use the central MARTY numerical policy.\n";
        outputFile << "\tauto hyperiso_regprop_it = param.realParams.find(\"reg_prop\");\n";
        outputFile << "\tbool hyperiso_has_regprop = hyperiso_regprop_it != param.realParams.end() && hyperiso_regprop_it->second != nullptr;\n";
        const bool read_split_sm_components = should_read_marty_split_sm_components(this->wilson);
        const bool expose_linker_components = (this->wilson == "CP10");
        outputFile << "\tdouble hyperiso_regprop_saved = " << default_regprop << ";\n";
        outputFile << "\tif (hyperiso_has_regprop) {\n";
        outputFile << "\t\thyperiso_regprop_saved = static_cast<double>(*hyperiso_regprop_it->second);\n";
        outputFile << "\t}\n";

        outputFile << "\tif (hyperiso_has_regprop) {\n";
        outputFile << "\t\t*hyperiso_regprop_it->second = " << default_regprop << ";\n";
        outputFile << "\t}\n";
        if (expose_linker_components) {
            // Use the full non-photon branch for the physical value.  The
            // VECTOR / SCALAR calls are diagnostics: they help identify where
            // CP10 comes from without risking double counting or losing boxes.
            outputFile << "\tauto hyperiso_bsm_non_photon = " + wilson + "(param);\n";
            outputFile << "\tauto hyperiso_bsm_vector = " + wilson + "_VECTOR(param);\n";
            outputFile << "\tauto hyperiso_bsm_scalar = " + wilson + "_SCALAR(param);\n";
        } else {
            outputFile << "\tauto hyperiso_bsm_non_photon = " + wilson + "(param);\n";
            outputFile << "\tauto hyperiso_bsm_vector = hyperiso_bsm_non_photon;\n";
            outputFile << "\tauto hyperiso_bsm_scalar = 0.0 * hyperiso_bsm_non_photon;\n";
        }
        if (read_split_sm_components) {
            if (expose_linker_components) {
                outputFile << "\tauto hyperiso_sm_non_photon = " + wilson + "_SM(param);\n";
                outputFile << "\tauto hyperiso_sm_vector = " + wilson + "_SM_VECTOR(param);\n";
                outputFile << "\tauto hyperiso_sm_scalar = " + wilson + "_SM_SCALAR(param);\n";
            } else {
                outputFile << "\tauto hyperiso_sm_non_photon = " + wilson + "_SM(param);\n";
                outputFile << "\tauto hyperiso_sm_vector = hyperiso_sm_non_photon;\n";
                outputFile << "\tauto hyperiso_sm_scalar = 0.0 * hyperiso_sm_non_photon;\n";
            }
        } else {
            outputFile << "\tauto hyperiso_sm_non_photon = 0.0 * hyperiso_bsm_non_photon;\n";
            outputFile << "\tauto hyperiso_sm_vector = 0.0 * hyperiso_bsm_vector;\n";
            outputFile << "\tauto hyperiso_sm_scalar = 0.0 * hyperiso_bsm_scalar;\n";
        }

        outputFile << "\tif (hyperiso_has_regprop) {\n";
        outputFile << "\t\t*hyperiso_regprop_it->second = " << photon_regprop << ";\n";
        outputFile << "\t}\n";
        outputFile << "\tauto hyperiso_bsm_photon = " + wilson + "_A(param);\n";
        if (read_split_sm_components) {
            outputFile << "\tauto hyperiso_sm_photon = " + wilson + "_SM_A(param);\n";
        } else {
            outputFile << "\tauto hyperiso_sm_photon = 0.0 * hyperiso_bsm_photon;\n";
        }

        outputFile << "\tif (hyperiso_has_regprop) {\n";
        outputFile << "\t\t*hyperiso_regprop_it->second = hyperiso_regprop_saved;\n";
        outputFile << "\t}\n";
        outputFile << "\tauto hyperiso_bsm_split = hyperiso_bsm_non_photon + hyperiso_bsm_photon;\n";
        outputFile << "\tauto hyperiso_sm_split = hyperiso_sm_non_photon + hyperiso_sm_photon;\n";
        outputFile << "\tauto hyperiso_total_split = hyperiso_sm_split + hyperiso_bsm_split;\n";
        if (full_target_generation) {
            outputFile << "\t// The analytical split contains the complete target model.\n";
            outputFile << "\twriteWilsonCoefficients(\"" + wilson + "\", hyperiso_bsm_split, Q_match, path);\n";
            outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_TARGET_SPLIT\", hyperiso_bsm_split, Q_match, path);\n";
        } else {
            outputFile << "\twriteWilsonCoefficients(\"" + wilson + "\", hyperiso_bsm_split, Q_match, path);\n";
        }
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_BSM_SPLIT\", hyperiso_bsm_split, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_SM_SPLIT\", hyperiso_sm_split, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_TOTAL_SPLIT\", hyperiso_total_split, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_NONPHOTON\", hyperiso_bsm_non_photon, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_A\", hyperiso_bsm_photon, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_SM_NONPHOTON\", hyperiso_sm_non_photon, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_SM_A\", hyperiso_sm_photon, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_VECTOR\", hyperiso_bsm_vector, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_SCALAR\", hyperiso_bsm_scalar, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_SM_VECTOR\", hyperiso_sm_vector, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_SM_SCALAR\", hyperiso_sm_scalar, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_TOTAL_VECTOR\", hyperiso_sm_vector + hyperiso_bsm_vector, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_TOTAL_SCALAR\", hyperiso_sm_scalar + hyperiso_bsm_scalar, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_SM_COMPONENT\", hyperiso_sm_split, Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_TOTAL_COMPONENT\", hyperiso_total_split, Q_match, path);\n";
    } else {
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "\", " + wilson + "(param), Q_match, path);\n";
    }

}

void MartyFileWriter::add_argpars(std::ofstream& outputFile) {
    const auto files = FileNameManager::getInstance(this->wilson, this->model);
    outputFile << "\tdouble Q_match = 80.379;\n";
    outputFile << "\tstd::string param_file_path = \"" << files->getParamFileName() << "\";\n";
    outputFile << "\tstd::string output_file_path = \"" << files->getCsvWilsonFileName() << "\";\n";
    outputFile << "\tbool raw_photon_diagnostic = false;\n";
    outputFile << "\tfor (int i = 1; i < argc; i++) {\n";
    outputFile << "\t\tif (std::string(argv[i]) == \"--Q_match\" || std::string(argv[i]) == \"-Q\") {\n";
    outputFile << "\t\t\tif (i + 1 >= argc) { throw std::runtime_error(\"Missing value after --Q_match\"); }\n";
    outputFile << "\t\t\tQ_match = std::stod(argv[i + 1]);\n";
    outputFile << "\t\t\ti++;\n";
    outputFile << "\t\t} else if (std::string(argv[i]) == \"--param-file\") {\n";
    outputFile << "\t\t\tif (i + 1 >= argc) { throw std::runtime_error(\"Missing value after --param-file\"); }\n";
    outputFile << "\t\t\tparam_file_path = argv[++i];\n";
    outputFile << "\t\t} else if (std::string(argv[i]) == \"--output-file\") {\n";
    outputFile << "\t\t\tif (i + 1 >= argc) { throw std::runtime_error(\"Missing value after --output-file\"); }\n";
    outputFile << "\t\t\toutput_file_path = argv[++i];\n";
    outputFile << "\t\t} else if (std::string(argv[i]) == \"--raw-photon-diagnostic\") {\n";
    outputFile << "\t\t\traw_photon_diagnostic = true;\n";
    outputFile << "\t\t} else if (std::string(argv[i]) == \"--help\" || std::string(argv[i]) == \"-h\") {\n";
    outputFile << "\t\t\tstd::cout << \"Options availables :\" << std::endl;\n";
    outputFile << "\t\t\tstd::cout << \"--Q_match/-Q : Value of Q_match (default 80.379)\" << std::endl;\n";
    outputFile << "\t\t\tstd::cout << \"--param-file : Per-invocation parameter CSV\" << std::endl;\n";
    outputFile << "\t\t\tstd::cout << \"--output-file : Per-invocation Wilson CSV\" << std::endl;\n";
    outputFile << "\t\t\tstd::cout << \"--raw-photon-diagnostic : Evaluate the unstable raw C9/CP9 photon-linker diagnostic\" << std::endl;\n";
    outputFile << "\t\t\tstd::cout << \"--help/-h : Affiche ce message.\" << std::endl;\n";
    outputFile << "\t\t\treturn 0;\n";
    outputFile << "\t\t}\n";
    outputFile << "\t}\n\n";
    
}

void MartyFileWriter::add_input_reader(std::ofstream& outputFile) {
    outputFile << "\tparam_t param;\n";
    outputFile << "\tstd::ifstream ParamFile(param_file_path);" << "\n";
    outputFile << "\tif (!ParamFile) { throw std::runtime_error(\"Cannot open MARTY parameter file: \" + param_file_path); }\n";
    outputFile << "\treadParams(ParamFile, param.realParams, param.complexParams);" << "\n";


}

#include "MartyInterface.h"
#include "ModelAPI.h"
#include "MartyParameterProxy.h"
#include "DefaultInterpreterPortsFactory.h"
#include "MartyRuntimeConfig.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <unordered_set>

namespace fs = std::filesystem;


MartyInterface::MartyInterface() {
    core_api = std::make_shared<ModelAPI>();
    param_proxy_sm = std::make_shared<MartyParameterProxy>(ParameterType::SM);
    param_proxy_bsm = std::make_shared<MartyParameterProxy>(ParameterType::BSM);
    ports = std::make_shared<DefaultInterpreterPortsFactory>();
}


void MartyInterface::compile_run(std::string wilson, std::string model) {
    if (!MartyRuntimeConfig::require_available("MartyInterface::compile_run").valid) {
        return;
    }

    GppCompilerStrategy compiler(model, wilson);
    if (!this->already_run(FileNameManager::getInstance(wilson, model)->getNumGeneratedFileName())){
        compiler.compile_run(FileNameManager::getInstance(wilson, model)->getGeneratedFileName(), FileNameManager::getInstance(wilson, model)->getExecutableFileName());
    }
}

void MartyInterface::generate(std::string wilson, std::string model, std::string model_path) {
    generate(std::move(wilson), model, model, std::move(model_path), false);
}

void MartyInterface::generate(std::string wilson,
                              std::string output_model,
                              std::string target_model,
                              std::string model_path,
                              bool sm_like_filter) {
    if (!MartyRuntimeConfig::require_available("MartyInterface::generate").valid) {
        return;
    }

    const auto model_template_index = resolve_model_template_index(target_model);
    invalidate_template_model_cache_if_needed(
        wilson, output_model, target_model, model_path, model_template_index, sm_like_filter
    );

    std::unique_ptr<ModelModifier> smModifier;
    smModifier = std::make_unique<GeneralModelModifier>(
        wilson, output_model, target_model, model_path, model_template_index, sm_like_filter
    );

    std::unique_ptr<TemplateManagerBase> templateManager = std::make_unique<NonNumericTemplateManager>(FileNameManager::getInstance(wilson, output_model)->getTemplateDir());
    templateManager->setModelAndWilson(output_model, wilson);
    templateManager->setModelModifier(std::move(smModifier));

    CodeGenerator codeGenerator(std::move(templateManager));

    codeGenerator.generate(wilson, FileNameManager::getInstance(wilson, output_model)->getGeneratedFileName());
}

void MartyInterface::generate_numlib(std::string wilson, std::string model) {
    generate_numlib(std::move(wilson), model, model);
}

void MartyInterface::generate_numlib(std::string wilson,
                                     std::string output_model,
                                     std::string target_model) {
    if (!MartyRuntimeConfig::require_available("MartyInterface::generate_numlib").valid) {
        return;
    }

    bool forceMode = false;
    auto file_names = FileNameManager::getInstance(wilson, output_model);
    const std::string cinematic_template = file_names->getGeneratedFileName();

    std::unique_ptr<SMParamSetter> sm_p_setter = std::make_unique<SMParamSetter>(
        target_model,
        specials_block,
        param_proxy_sm,
        param_proxy_bsm,
        cinematic_template
    );

    std::unique_ptr<GeneralNumModelModifier> ModelModifier = std::make_unique<GeneralNumModelModifier>(
        wilson,
        output_model,
        target_model,
        std::move(sm_p_setter),
        core_api,
        ports,
        forceMode
    );
    
    std::unique_ptr<TemplateManagerBase> templateManager = std::make_unique<NumericTemplateManager>(file_names->getLibDir());
    templateManager->setModelAndWilson(output_model, wilson);
    templateManager->setNumModelModifier(std::move(ModelModifier));
    this->dependencies.insert_or_assign(wilson, templateManager->get_dependencies());
    CodeGenerator codeGenerator(std::move(templateManager));
    std::string file_path = file_names->getNumGeneratedFileName();
    codeGenerator.generate(file_path, file_path);
}

void MartyInterface::compile_run_libs(std::string wilson, std::string model, double Q_match) {
    if (!MartyRuntimeConfig::require_available("MartyInterface::compile_run_libs").valid) {
        return;
    }

    MakeCompilerStrategy compiler(model, wilson);
    compiler.set_Q_match(Q_match);
    compiler.compile_run(FileNameManager::getInstance(wilson, model)->getLibDir(), FileNameManager::getInstance(wilson,model)->getNumExecutableFileName());
}

void MartyInterface::calculate(std::string wilson, std::string model, double Q_match, std::string model_path) {
    calculate(std::move(wilson), model, model, Q_match, std::move(model_path), false);
}

void MartyInterface::calculate(std::string wilson,
                               std::string output_model,
                               std::string target_model,
                               double Q_match,
                               std::string model_path,
                               bool sm_like_filter) {
    if (!MartyRuntimeConfig::require_available("MartyInterface::calculate").valid) {
        return;
    }

    generate(wilson, output_model, target_model, model_path, sm_like_filter);
    compile_run(wilson, output_model);
    generate_numlib(wilson, output_model, target_model);
    compile_run_libs(wilson, output_model, Q_match);
}


std::optional<int> MartyInterface::resolve_model_template_index(const std::string& model) const {
    std::string model_upper = model;
    std::transform(model_upper.begin(), model_upper.end(), model_upper.begin(), [](unsigned char c) {
        return static_cast<char>(std::toupper(c));
    });

    if (model_upper != "THDM") {
        return std::nullopt;
    }

    if (!param_proxy_bsm) {
        LOG_ERROR("MartyConfigError", "Cannot instantiate the templated THDM MARTY model: no BSM parameter proxy is available to read MINPAR(24). ",
                  "Set the THDM Yukawa type in the LHA card or provide a BSM parameter provider before MARTY generation.");
    }

    const double raw_type = (*param_proxy_bsm)("MINPAR", LhaID(24));
    const int type = static_cast<int>(std::lround(raw_type));

    if (std::abs(raw_type - static_cast<double>(type)) > 1e-9 || type < 1 || type > 4) {
        LOG_ERROR("MartyConfigError", "Invalid THDM Yukawa type MINPAR(24)=", raw_type,
                  ". MARTY THDM generation expects an integer type in {1,2,3,4}.");
    }

    LOG_INFO("MartyInterface", "Using THDM Yukawa type ", type, " from MINPAR(24) for MARTY generation.");
    return type;
}

void MartyInterface::invalidate_template_model_cache_if_needed(const std::string& wilson,
                                                               const std::string& output_model,
                                                               const std::string& target_model,
                                                               const std::string& model_path,
                                                               std::optional<int> model_template_index,
                                                               bool sm_like_filter) const {
    const auto files = FileNameManager::getInstance(wilson, output_model);
    const std::string expected_signature = GeneralModelModifier::modelSignature(
        target_model,
        model_path,
        model_template_index
    );

    const std::unordered_set<std::string> semileptonic_templates = {"C9", "C10", "CP9", "CP10"};
    const bool needs_semileptonic_template_abi = semileptonic_templates.contains(wilson);
    constexpr const char* expected_semileptonic_template_abi =
        "HYPERISO_MARTY_TEMPLATE_ABI: semileptonic-local-vertex-diagnostics-v1";

    bool stale = true;
    {
        std::ifstream in(files->getGeneratedFileName());
        std::string line;
        bool has_expected_signature = false;
        bool has_expected_filter_mode = !sm_like_filter;
        bool has_expected_template_abi = !needs_semileptonic_template_abi;
        while (std::getline(in, line)) {
            if (line.find(expected_signature) != std::string::npos) {
                has_expected_signature = true;
            }
            if (line.find("HYPERISO_MARTY_SM_LIKE_FILTER:") != std::string::npos) {
                has_expected_filter_mode = sm_like_filter;
            }
            if (needs_semileptonic_template_abi
                && line.find(expected_semileptonic_template_abi) != std::string::npos) {
                has_expected_template_abi = true;
            }
            if (has_expected_signature && has_expected_filter_mode && has_expected_template_abi) {
                stale = false;
                break;
            }
            if (line.find("HYPERISO_MARTY_MODEL_SIGNATURE:") != std::string::npos
                && !has_expected_signature) {
                break;
            }
        }
    }

    if (!stale) {
        return;
    }

    std::error_code ec;
    fs::remove(files->getGeneratedFileName(), ec);
    ec.clear();
    fs::remove(files->getExecutableFileName(), ec);
    ec.clear();
    fs::remove_all(files->getLibDir(), ec);
    ec.clear();
    fs::remove(files->getCsvWilsonFileName(), ec);

    LOG_INFO("MartyInterface", "Invalidated stale MARTY cache for ", wilson, " / ", output_model,
             " because the model signature is now ", expected_signature);
}

std::unordered_set<InterpretedParam> MartyInterface::get_dependencies(std::string wilson) {
    if (!this->dependencies.contains(wilson)) {
        LOG_ERROR("KeyError", "Trying to access dependencies for unknown wilson coefficient", wilson, "in WilsonInterface.");
    }
    
    return this->dependencies.at(wilson);
}

bool MartyInterface::already_run(std::string&& outputBinary) {
    struct stat buffer;
    if (stat(outputBinary.c_str(), &buffer) != 0) {
        return false;
    }
    if (buffer.st_size == 0) {
        return false;
    }
    LOG_DEBUG("Already run !");
    return true;
}

std::string MartyInterface::output_binary_name(std::string& wilson, std::string& model) {
        return "generated_" + wilson+"_" + model + ".cpp";
    }

std::set<std::string> MartyInterface::get_special_blocks() {
    return this->specials_block;
}
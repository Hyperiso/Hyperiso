#include "MartyInterface.h"
#include "ModelAPI.h"
#include "MartyParameterProxy.h"
#include "DefaultInterpreterPortsFactory.h"
#include "MartyRuntimeConfig.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>

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
    if (!MartyRuntimeConfig::require_available("MartyInterface::generate").valid) {
        return;
    }

    const auto model_template_index = resolve_model_template_index(model);
    invalidate_template_model_cache_if_needed(wilson, model, model_path, model_template_index);

    std::unique_ptr<ModelModifier> smModifier;
    smModifier = std::make_unique<GeneralModelModifier>(wilson, model, model_path, model_template_index);

    std::unique_ptr<TemplateManagerBase> templateManager = std::make_unique<NonNumericTemplateManager>(FileNameManager::getInstance(wilson, model)->getTemplateDir());
    templateManager->setModelAndWilson(model, wilson);
    templateManager->setModelModifier(std::move(smModifier));

    CodeGenerator codeGenerator(std::move(templateManager));

    codeGenerator.generate(wilson, FileNameManager::getInstance(wilson, model)->getGeneratedFileName());
}

void MartyInterface::generate_numlib(std::string wilson, std::string model) {
    if (!MartyRuntimeConfig::require_available("MartyInterface::generate_numlib").valid) {
        return;
    }

    bool forceMode = false;
    auto file_names = FileNameManager::getInstance(wilson, model);
    const std::string cinematic_template = file_names->getGeneratedFileName();

    std::unique_ptr<SMParamSetter> sm_p_setter = std::make_unique<SMParamSetter>(
        model,
        specials_block,
        param_proxy_sm,
        param_proxy_bsm,
        cinematic_template
    );
    std::unique_ptr<GeneralNumModelModifier> ModelModifier = std::make_unique<GeneralNumModelModifier>(wilson, model, std::move(sm_p_setter), core_api, ports, forceMode);
    
    std::unique_ptr<TemplateManagerBase> templateManager = std::make_unique<NumericTemplateManager>(file_names->getLibDir());
    templateManager->setModelAndWilson(model, wilson);
    templateManager->setNumModelModifier(std::move(ModelModifier));
    this->dependencies.emplace(wilson, templateManager->get_dependencies());
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
    if (!MartyRuntimeConfig::require_available("MartyInterface::calculate").valid) {
        return;
    }

    generate(wilson, model, model_path);
    compile_run(wilson, model);
    generate_numlib(wilson, model);
    compile_run_libs(wilson, model, Q_match);
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
                                                               const std::string& model,
                                                               const std::string& model_path,
                                                               std::optional<int> model_template_index) const {
    const auto files = FileNameManager::getInstance(wilson, model);
    const std::string expected_signature = GeneralModelModifier::modelSignature(
        model,
        model_path,
        model_template_index
    );

    bool stale = true;
    {
        std::ifstream in(files->getGeneratedFileName());
        std::string line;
        while (std::getline(in, line)) {
            if (line.find(expected_signature) != std::string::npos) {
                stale = false;
                break;
            }
            if (line.find("HYPERISO_MARTY_MODEL_SIGNATURE:") != std::string::npos) {
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

    LOG_INFO("MartyInterface", "Invalidated stale MARTY cache for ", wilson, " / ", model,
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
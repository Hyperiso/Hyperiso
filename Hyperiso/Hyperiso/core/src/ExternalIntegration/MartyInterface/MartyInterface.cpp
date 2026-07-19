#include "MartyInterface.h"
#include "ModelAPI.h"
#include "MartyParameterProxy.h"
#include "DefaultInterpreterPortsFactory.h"
#include "MartyRuntimeConfig.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>

namespace fs = std::filesystem;

namespace {
std::string template_signature(const std::string& wilson,
                               const std::shared_ptr<FileNameManager>& files);
std::string generation_mode_marker(const std::string& wilson,
                                   bool sm_like_filter,
                                   bool bsm_only_generation,
                                   bool full_target_generation);
void append_cache_metadata_if_missing(const fs::path& generated_file,
                                      const std::string& model_signature,
                                      const std::string& template_signature_value,
                                      const std::string& mode_marker);
} // namespace


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
    const auto files = FileNameManager::getInstance(wilson, model);
    if (!this->already_run(files->getExecutableFileName())) {
        compiler.compile_run(files->getGeneratedFileName(), files->getExecutableFileName());
    }
}

void MartyInterface::generate(std::string wilson, std::string model, std::string model_path) {
    generate(std::move(wilson), model, model, std::move(model_path), false, false, false);
}

void MartyInterface::generate(std::string wilson,
                              std::string output_model,
                              std::string target_model,
                              std::string model_path,
                              bool sm_like_filter,
                              bool bsm_split_generation,
                              bool full_target_generation) {
    if (!MartyRuntimeConfig::require_available("MartyInterface::generate").valid) {
        return;
    }

    const auto model_template_index = resolve_model_template_index(target_model);
    invalidate_template_model_cache_if_needed(
        wilson, output_model, target_model, model_path, model_template_index,
        sm_like_filter, bsm_split_generation, full_target_generation
    );

    std::unique_ptr<ModelModifier> smModifier;
    smModifier = std::make_unique<GeneralModelModifier>(
        wilson, output_model, target_model, model_path, model_template_index,
        sm_like_filter, bsm_split_generation, full_target_generation
    );

    const auto files = FileNameManager::getInstance(wilson, output_model);
    std::unique_ptr<TemplateManagerBase> templateManager = std::make_unique<NonNumericTemplateManager>(files->getTemplateDir());
    templateManager->setModelAndWilson(output_model, wilson);
    templateManager->setModelModifier(std::move(smModifier));

    CodeGenerator codeGenerator(std::move(templateManager));

    codeGenerator.generate(wilson, files->getGeneratedFileName());
    append_cache_metadata_if_missing(
        files->getGeneratedFileName(),
        GeneralModelModifier::modelSignature(target_model, model_path, model_template_index),
        template_signature(wilson, files),
        generation_mode_marker(wilson, sm_like_filter, bsm_split_generation, full_target_generation)
    );
}

void MartyInterface::generate_numlib(std::string wilson, std::string model) {
    generate_numlib(std::move(wilson), model, model, false, false);
}

void MartyInterface::generate_numlib(std::string wilson,
                                     std::string output_model,
                                     std::string target_model,
                                     bool bsm_split_generation,
                                     bool full_target_generation) {
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
        forceMode,
        bsm_split_generation,
        full_target_generation
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
    calculate(std::move(wilson), model, model, Q_match, std::move(model_path), false, false, false);
}

void MartyInterface::calculate(std::string wilson,
                               std::string output_model,
                               std::string target_model,
                               double Q_match,
                               std::string model_path,
                               bool sm_like_filter,
                               bool bsm_split_generation,
                               bool full_target_generation) {
    if (!MartyRuntimeConfig::require_available("MartyInterface::calculate").valid) {
        return;
    }

    generate(
        wilson, output_model, target_model, model_path,
        sm_like_filter, bsm_split_generation, full_target_generation
    );
    compile_run(wilson, output_model);
    generate_numlib(
        wilson, output_model, target_model,
        bsm_split_generation, full_target_generation
    );
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

namespace {

constexpr const char* kMartyCacheAbi = "HYPERISO_MARTY_CACHE_ABI: pyhyperiso-1.0.3-v2";

std::string stable_file_fingerprint(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("Cannot fingerprint MARTY template file: " + path.string());
    }

    std::uint64_t hash = 14695981039346656037ULL;
    char buffer[8192];
    while (input.read(buffer, sizeof(buffer)) || input.gcount() > 0) {
        const auto count = input.gcount();
        for (std::streamsize i = 0; i < count; ++i) {
            hash ^= static_cast<unsigned char>(buffer[i]);
            hash *= 1099511628211ULL;
        }
    }

    std::ostringstream result;
    result << std::hex << std::setw(16) << std::setfill('0') << hash;
    return result.str();
}

std::string normalized_path(const fs::path& path) {
    std::error_code ec;
    fs::path normalized = fs::weakly_canonical(path, ec);
    if (ec) {
        ec.clear();
        normalized = fs::absolute(path, ec);
    }
    return normalized.lexically_normal().string();
}

std::string template_signature(const std::string& wilson,
                               const std::shared_ptr<FileNameManager>& files) {
    const fs::path path = fs::path(files->getTemplateDir()) / (wilson + ".cpp");
    return "HYPERISO_MARTY_TEMPLATE_SIGNATURE: path=" + normalized_path(path)
         + "; fnv1a64=" + stable_file_fingerprint(path);
}

bool uses_split_regprop_policy(const std::string& wilson) {
    return wilson == "C9" || wilson == "CP9" || wilson == "CP10";
}

std::string generation_mode(const std::string& wilson,
                            bool sm_like_filter,
                            bool bsm_only_generation,
                            bool full_target_generation) {
    if (sm_like_filter) {
        return "sm-like";
    }
    if (full_target_generation && bsm_only_generation && uses_split_regprop_policy(wilson)) {
        return "target-regprop-split";
    }
    if (full_target_generation) {
        return "target-full";
    }
    if (bsm_only_generation && uses_split_regprop_policy(wilson)) {
        return "bsm-regprop-split";
    }
    if (bsm_only_generation) {
        return "bsm-only";
    }
    return "full";
}

std::string generation_mode_marker(const std::string& wilson,
                                   bool sm_like_filter,
                                   bool bsm_only_generation,
                                   bool full_target_generation) {
    return "HYPERISO_MARTY_GENERATION_MODE: "
         + generation_mode(
             wilson,
             sm_like_filter,
             bsm_only_generation,
             full_target_generation
         );
}

void append_cache_metadata_if_missing(const fs::path& generated_file,
                                      const std::string& model_signature,
                                      const std::string& template_signature_value,
                                      const std::string& mode_marker) {
    std::ifstream input(generated_file);
    if (!input) {
        throw std::runtime_error(
            "MARTY source generation did not create the expected file: " + generated_file.string()
        );
    }

    bool has_cache_abi = false;
    std::string line;
    while (std::getline(input, line)) {
        if (line.find(kMartyCacheAbi) != std::string::npos) {
            has_cache_abi = true;
            break;
        }
    }
    if (has_cache_abi) {
        return;
    }

    std::ofstream output(generated_file, std::ios::app);
    if (!output) {
        throw std::runtime_error(
            "Cannot append MARTY cache metadata to: " + generated_file.string()
        );
    }
    output << "\n// " << kMartyCacheAbi << "\n";
    output << "// " << model_signature << "\n";
    output << "// " << template_signature_value << "\n";
    output << "// " << mode_marker << "\n";
}

} // namespace

void MartyInterface::invalidate_template_model_cache_if_needed(const std::string& wilson,
                                                               const std::string& output_model,
                                                               const std::string& target_model,
                                                               const std::string& model_path,
                                                               std::optional<int> model_template_index,
                                                               bool sm_like_filter,
                                                               bool bsm_split_generation,
                                                               bool full_target_generation) const {
    const auto files = FileNameManager::getInstance(wilson, output_model);
    const std::string expected_model_signature = GeneralModelModifier::modelSignature(
        target_model,
        model_path,
        model_template_index
    );
    const std::string expected_template_signature = template_signature(wilson, files);
    const std::string expected_mode = generation_mode_marker(
        wilson,
        sm_like_filter,
        bsm_split_generation,
        full_target_generation
    );

    bool file_present = false;
    bool has_cache_abi = false;
    bool has_model_signature = false;
    bool has_template_signature = false;
    bool has_generation_mode = false;

    {
        std::ifstream in(files->getGeneratedFileName());
        file_present = static_cast<bool>(in);
        std::string line;
        while (std::getline(in, line)) {
            has_cache_abi = has_cache_abi || line.find(kMartyCacheAbi) != std::string::npos;
            has_model_signature = has_model_signature || line.find(expected_model_signature) != std::string::npos;
            has_template_signature = has_template_signature || line.find(expected_template_signature) != std::string::npos;
            has_generation_mode = has_generation_mode || line.find(expected_mode) != std::string::npos;
        }
    }

    const bool stale = !file_present
                    || !has_cache_abi
                    || !has_model_signature
                    || !has_template_signature
                    || !has_generation_mode;

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

    std::string reason;
    if (!file_present) {
        reason = "generated file is missing";
    } else if (!has_cache_abi) {
        reason = "cache ABI mismatch";
    } else if (!has_model_signature) {
        reason = "model path/content signature mismatch";
    } else if (!has_template_signature) {
        reason = "template content signature mismatch";
    } else if (!has_generation_mode) {
        reason = "generation mode mismatch";
    } else {
        reason = "cache metadata mismatch";
    }

    LOG_INFO("MartyInterface", "Invalidated stale MARTY cache for ", wilson, " / ", output_model,
             " because ", reason, ". Expected mode: ", expected_mode,
             "; expected model signature: ", expected_model_signature);
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

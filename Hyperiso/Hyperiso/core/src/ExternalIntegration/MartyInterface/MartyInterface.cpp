#include "MartyInterface.h"
#include "ModelAPI.h"
#include "MartyParameterProxy.h"
#include "DefaultInterpreterPortsFactory.h"
#include "MartyRuntimeConfig.h"
#include "MartyAdapter.h"
#include "ParamWriter.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <map>
#include <mutex>
#include <optional>
#include <random>
#include <regex>
#include <shared_mutex>
#include <sstream>
#include <thread>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

namespace {
std::shared_mutex marty_artifact_mutex;
std::mutex marty_legacy_csv_mutex;
std::atomic<std::uint64_t> marty_run_counter {0};

constexpr const char* kMartyCacheAbi = "HYPERISO_MARTY_CACHE_ABI: pyhyperiso-1.0.4-v1";

std::string sanitize_path_component(std::string value) {
    for (char& c : value) {
        const auto uc = static_cast<unsigned char>(c);
        if (!std::isalnum(uc) && c != '-' && c != '_') {
            c = '_';
        }
    }
    return value.empty() ? std::string("unnamed") : value;
}

fs::path make_invocation_directory(const std::shared_ptr<FileNameManager>& files,
                                   const std::string& wilson,
                                   const std::string& model) {
    static const std::uint64_t process_nonce = [] {
        std::random_device random;
        return (static_cast<std::uint64_t>(random()) << 32)
             ^ static_cast<std::uint64_t>(random());
    }();
    const auto counter = marty_run_counter.fetch_add(1, std::memory_order_relaxed);
    const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
    const auto thread_hash = std::hash<std::thread::id>{}(std::this_thread::get_id());

    std::ostringstream name;
    name << sanitize_path_component(model) << "_"
         << sanitize_path_component(wilson) << "_"
         << std::hex << process_nonce << "_" << now << "_"
         << thread_hash << "_" << counter;

    const fs::path dir = fs::path(files->getOutputDir()) / "runs" / name.str();
    std::error_code ec;
    fs::create_directories(dir, ec);
    if (ec) {
        throw std::runtime_error(
            "Cannot create invocation-local MARTY directory: " + dir.string()
            + " (" + ec.message() + ")"
        );
    }
    return dir;
}

void write_parameter_snapshot(const fs::path& path,
                              const std::unordered_map<std::string, double>& params) {
    std::ofstream output(path, std::ios::trunc);
    if (!output) {
        throw std::runtime_error("Cannot write MARTY parameter snapshot: " + path.string());
    }

    // Reuse the exact legacy serializer instead of changing the numerical
    // representation merely because the file is invocation-local.  The path
    // isolation provides thread safety; parameter rounding must remain
    // backward-compatible with the pre-thread-safe implementation.
    ParamWriter parameter_writer;
    parameter_writer.writeParams(output, params);
    output.flush();
    if (!output) {
        throw std::runtime_error("Failed while writing MARTY parameter snapshot: " + path.string());
    }
}

void publish_legacy_csv(const fs::path& isolated, const fs::path& legacy) {
    std::lock_guard<std::mutex> lock(marty_legacy_csv_mutex);

    const auto split_csv_line = [](const std::string& line) {
        std::vector<std::string> cells;
        std::stringstream stream(line);
        std::string cell;
        while (std::getline(stream, cell, ',')) {
            cells.push_back(cell);
        }
        return cells;
    };

    struct CsvTable {
        std::vector<std::string> headers;
        std::vector<std::vector<std::string>> rows;
    };

    const auto read_table = [&](const fs::path& path, bool required) {
        CsvTable table;
        std::ifstream input(path);
        if (!input) {
            if (required) {
                throw std::runtime_error("Cannot read MARTY CSV: " + path.string());
            }
            return table;
        }
        std::string line;
        if (std::getline(input, line)) {
            table.headers = split_csv_line(line);
        }
        while (std::getline(input, line)) {
            if (!line.empty()) {
                table.rows.push_back(split_csv_line(line));
            }
        }
        return table;
    };

    CsvTable incoming = read_table(isolated, true);
    if (incoming.headers.empty() || incoming.headers.front() != "Q_match") {
        throw std::runtime_error("Invalid invocation-local MARTY CSV: " + isolated.string());
    }

    CsvTable merged = read_table(legacy, false);
    if (merged.headers.empty()) {
        merged.headers = {"Q_match"};
    }
    if (merged.headers.front() != "Q_match") {
        throw std::runtime_error("Invalid legacy MARTY CSV: " + legacy.string());
    }

    std::unordered_map<std::string, std::size_t> merged_columns;
    for (std::size_t i = 0; i < merged.headers.size(); ++i) {
        merged_columns.emplace(merged.headers[i], i);
    }
    for (std::size_t i = 1; i < incoming.headers.size(); ++i) {
        if (!merged_columns.contains(incoming.headers[i])) {
            merged_columns.emplace(incoming.headers[i], merged.headers.size());
            merged.headers.push_back(incoming.headers[i]);
            for (auto& row : merged.rows) {
                row.resize(merged.headers.size(), "NaN");
            }
        }
    }

    std::map<double, std::size_t> merged_rows;
    for (std::size_t i = 0; i < merged.rows.size(); ++i) {
        merged.rows[i].resize(merged.headers.size(), "NaN");
        if (!merged.rows[i].empty()) {
            merged_rows[std::stod(merged.rows[i][0])] = i;
        }
    }

    for (const auto& incoming_row : incoming.rows) {
        if (incoming_row.empty()) {
            continue;
        }
        const double q_match = std::stod(incoming_row[0]);
        std::size_t row_index = 0;
        const auto existing = merged_rows.find(q_match);
        if (existing == merged_rows.end()) {
            row_index = merged.rows.size();
            merged.rows.emplace_back(merged.headers.size(), "NaN");
            merged.rows.back()[0] = incoming_row[0];
            merged_rows.emplace(q_match, row_index);
        } else {
            row_index = existing->second;
        }

        for (std::size_t i = 1; i < incoming.headers.size() && i < incoming_row.size(); ++i) {
            merged.rows[row_index][merged_columns.at(incoming.headers[i])] = incoming_row[i];
        }
    }

    std::error_code ec;
    fs::create_directories(legacy.parent_path(), ec);
    if (ec) {
        throw std::runtime_error("Cannot create MARTY CSV directory: " + ec.message());
    }

    const fs::path tmp = legacy.string() + ".tmp." + std::to_string(
        marty_run_counter.fetch_add(1, std::memory_order_relaxed)
    );
    {
        std::ofstream output(tmp, std::ios::trunc);
        if (!output) {
            throw std::runtime_error("Cannot write temporary MARTY CSV: " + tmp.string());
        }
        for (std::size_t i = 0; i < merged.headers.size(); ++i) {
            output << merged.headers[i] << (i + 1 == merged.headers.size() ? '\n' : ',');
        }
        for (const auto& row : merged.rows) {
            for (std::size_t i = 0; i < merged.headers.size(); ++i) {
                output << (i < row.size() ? row[i] : "NaN")
                       << (i + 1 == merged.headers.size() ? '\n' : ',');
            }
        }
        output.flush();
        if (!output) {
            std::error_code cleanup_ec;
            fs::remove(tmp, cleanup_ec);
            throw std::runtime_error("Failed while writing temporary MARTY CSV: " + tmp.string());
        }
    }

    fs::rename(tmp, legacy, ec);
    if (ec) {
        std::error_code cleanup_ec;
        fs::remove(tmp, cleanup_ec);
        throw std::runtime_error("Cannot atomically publish MARTY CSV: " + ec.message());
    }
}

std::string template_signature(const std::string& wilson,
                               const std::shared_ptr<FileNameManager>& files);
MartyOrderPolicy effective_order_policy(bool sm_like_filter,
                                          bool bsm_only_generation,
                                          bool full_target_generation);
std::string generation_mode_marker(const std::string& wilson,
                                   bool sm_like_filter,
                                   bool bsm_only_generation,
                                   bool full_target_generation,
                                   MartyOrderPolicy order_policy);
void append_cache_metadata_if_missing(const fs::path& generated_file,
                                      const std::string& model_signature,
                                      const std::string& template_signature_value,
                                      const std::string& mode_marker);
bool template_needs_generic_tree_first(const std::string& wilson,
                                       const std::shared_ptr<FileNameManager>& files);
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
    const auto files = FileNameManager::getInstance(wilson, output_model);
    const bool tree_first_fallback = template_needs_generic_tree_first(wilson, files);
    const MartyOrderPolicy order_policy = effective_order_policy(
        sm_like_filter, bsm_split_generation, full_target_generation
    );
    invalidate_template_model_cache_if_needed(
        wilson, output_model, target_model, model_path, model_template_index,
        sm_like_filter, bsm_split_generation, full_target_generation
    );

    std::unique_ptr<ModelModifier> smModifier;
    smModifier = std::make_unique<GeneralModelModifier>(
        wilson, output_model, target_model, model_path, model_template_index,
        sm_like_filter, bsm_split_generation, full_target_generation,
        tree_first_fallback, order_policy
    );

    std::unique_ptr<TemplateManagerBase> templateManager = std::make_unique<NonNumericTemplateManager>(files->getTemplateDir());
    templateManager->setModelAndWilson(output_model, wilson);
    templateManager->setModelModifier(std::move(smModifier));

    CodeGenerator codeGenerator(std::move(templateManager));

    codeGenerator.generate(wilson, files->getGeneratedFileName());
    append_cache_metadata_if_missing(
        files->getGeneratedFileName(),
        GeneralModelModifier::modelSignature(target_model, model_path, model_template_index),
        template_signature(wilson, files),
        generation_mode_marker(
            wilson, sm_like_filter, bsm_split_generation, full_target_generation, order_policy
        )
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
    const auto discovered_dependencies = templateManager->get_dependencies();
    auto& cached_dependencies = this->dependencies[wilson];
    cached_dependencies.insert(discovered_dependencies.begin(), discovered_dependencies.end());
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

    const fs::path isolated_csv = calculate_isolated(
        wilson,
        output_model,
        target_model,
        Q_match,
        model_path,
        sm_like_filter,
        bsm_split_generation,
        full_target_generation
    );
    if (isolated_csv.empty()) {
        return;
    }

    const fs::path legacy_csv = FileNameManager::getInstance(wilson, output_model)->getCsvWilsonFileName();
    try {
        publish_legacy_csv(isolated_csv, legacy_csv);
    } catch (...) {
        std::error_code cleanup_ec;
        fs::remove_all(isolated_csv.parent_path(), cleanup_ec);
        throw;
    }

    std::error_code cleanup_ec;
    fs::remove_all(isolated_csv.parent_path(), cleanup_ec);
}

void MartyInterface::compile_numlib(const std::string& wilson, const std::string& model) {
    const auto files = FileNameManager::getInstance(wilson, model);
    MakeCompilerStrategy compiler(model, wilson);
    if (!compiler.check_if_compile(files->getNumExecutableFileName())) {
        compiler.compile(files->getLibDir(), files->getNumExecutableFileName());
    }
}

std::unordered_map<std::string, double> MartyInterface::snapshot_numeric_params(
    const std::string& wilson,
    const std::string& output_model,
    const std::string& target_model,
    bool bsm_split_generation,
    bool full_target_generation
) {
    const auto files = FileNameManager::getInstance(wilson, output_model);
    auto setter = std::make_unique<SMParamSetter>(
        target_model,
        specials_block,
        param_proxy_sm,
        param_proxy_bsm,
        files->getGeneratedFileName()
    );
    GeneralNumModelModifier modifier(
        wilson,
        output_model,
        target_model,
        std::move(setter),
        core_api,
        ports,
        false,
        bsm_split_generation,
        full_target_generation
    );
    return modifier.get_params();
}

void MartyInterface::ensure_built(const std::string& wilson,
                                  const std::string& output_model,
                                  const std::string& target_model,
                                  const std::string& model_path,
                                  bool sm_like_filter,
                                  bool bsm_split_generation,
                                  bool full_target_generation) {
    generate(
        wilson,
        output_model,
        target_model,
        model_path,
        sm_like_filter,
        bsm_split_generation,
        full_target_generation
    );
    compile_run(wilson, output_model);
    generate_numlib(
        wilson,
        output_model,
        target_model,
        bsm_split_generation,
        full_target_generation
    );
    compile_numlib(wilson, output_model);
}

bool MartyInterface::artifacts_ready(const std::string& wilson,
                                     const std::string& output_model,
                                     const std::string& target_model,
                                     const std::string& model_path,
                                     bool sm_like_filter,
                                     bool bsm_split_generation,
                                     bool full_target_generation) const {
    const auto files = FileNameManager::getInstance(wilson, output_model);
    const auto model_template_index = resolve_model_template_index(target_model);
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
        full_target_generation,
        effective_order_policy(sm_like_filter, bsm_split_generation, full_target_generation)
    );

    std::ifstream generated(files->getGeneratedFileName());
    if (!generated) {
        return false;
    }

    bool has_cache_abi = false;
    bool has_model_signature = false;
    bool has_template_signature = false;
    bool has_generation_mode = false;
    std::string line;
    while (std::getline(generated, line)) {
        has_cache_abi = has_cache_abi || line.find(kMartyCacheAbi) != std::string::npos;
        has_model_signature = has_model_signature || line.find(expected_model_signature) != std::string::npos;
        has_template_signature = has_template_signature || line.find(expected_template_signature) != std::string::npos;
        has_generation_mode = has_generation_mode || line.find(expected_mode) != std::string::npos;
    }

    const auto non_empty_file = [](const fs::path& path) {
        std::error_code ec;
        return fs::is_regular_file(path, ec) && !ec && fs::file_size(path, ec) > 0 && !ec;
    };
    const auto has_generation_marker = [](const fs::path& path) {
        std::ifstream input(path);
        std::string header;
        for (int i = 0; i < 16 && std::getline(input, header); ++i) {
            if (header.find("//42") != std::string::npos) {
                return true;
            }
        }
        return false;
    };

    return has_cache_abi
        && has_model_signature
        && has_template_signature
        && has_generation_mode
        && non_empty_file(files->getExecutableFileName())
        && non_empty_file(files->getNumGeneratedFileName())
        && has_generation_marker(files->getNumGeneratedFileName())
        && non_empty_file(files->getNumExecutableFileName())
        && this->dependencies.contains(wilson);
}

std::string MartyInterface::calculate_isolated(std::string wilson,
                                               std::string output_model,
                                               std::string target_model,
                                               double Q_match,
                                               std::string model_path,
                                               bool sm_like_filter,
                                               bool bsm_split_generation,
                                               bool full_target_generation) {
    if (!MartyRuntimeConfig::require_available("MartyInterface::calculate_isolated").valid) {
        return {};
    }

    const auto execute_isolated = [&]() -> std::string {
        const auto files = FileNameManager::getInstance(wilson, output_model);
        const fs::path run_dir = make_invocation_directory(files, wilson, output_model);
        const fs::path param_file = run_dir / "paramlist.csv";
        const fs::path output_file = run_dir / "wilson.csv";

        try {
            write_parameter_snapshot(
                param_file,
                snapshot_numeric_params(
                    wilson,
                    output_model,
                    target_model,
                    bsm_split_generation,
                    full_target_generation
                )
            );

            MakeCompilerStrategy compiler(output_model, wilson);
            compiler.set_Q_match(Q_match);
            compiler.set_param_file(param_file);
            compiler.set_output_file(output_file);
            compiler.compile_run(files->getLibDir(), files->getNumExecutableFileName());

            std::error_code ec;
            if (!fs::is_regular_file(output_file, ec) || ec || fs::file_size(output_file, ec) == 0 || ec) {
                throw std::runtime_error(
                    "MARTY numeric execution did not create a non-empty invocation-local CSV: "
                    + output_file.string()
                );
            }
            return output_file.string();
        } catch (...) {
            std::error_code cleanup_ec;
            fs::remove_all(run_dir, cleanup_ec);
            throw;
        }
    };

    {
        std::shared_lock<std::shared_mutex> read_lock(marty_artifact_mutex);
        if (artifacts_ready(
                wilson,
                output_model,
                target_model,
                model_path,
                sm_like_filter,
                bsm_split_generation,
                full_target_generation
            )) {
            return execute_isolated();
        }
    }

    std::unique_lock<std::shared_mutex> build_lock(marty_artifact_mutex);
    if (!artifacts_ready(
            wilson,
            output_model,
            target_model,
            model_path,
            sm_like_filter,
            bsm_split_generation,
            full_target_generation
        )) {
        ensure_built(
            wilson,
            output_model,
            target_model,
            model_path,
            sm_like_filter,
            bsm_split_generation,
            full_target_generation
        );
    }
    return execute_isolated();
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

bool template_needs_generic_tree_first(const std::string& wilson,
                                       const std::shared_ptr<FileNameManager>& files) {
    if (uses_split_regprop_policy(wilson)) {
        return false;
    }

    const fs::path path = fs::path(files->getTemplateDir()) / (wilson + ".cpp");
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Cannot inspect MARTY template order policy: " + path.string());
    }

    const std::string source(
        (std::istreambuf_iterator<char>(input)),
        std::istreambuf_iterator<char>()
    );
    static const std::regex tree_call(
        R"(computeWilsonCoefficients\s*\(\s*(?:mty::Order::)?TreeLevel)"
    );
    static const std::regex loop_call(
        R"(computeWilsonCoefficients\s*\(\s*(?:mty::Order::)?OneLoop)"
    );

    return std::regex_search(source, loop_call) && !std::regex_search(source, tree_call);
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

MartyOrderPolicy effective_order_policy(bool sm_like_filter,
                                          bool bsm_only_generation,
                                          bool full_target_generation) {
    // The user policy belongs to the configured BSM target.  The separately
    // generated SM baseline must retain AUTO so that tree-level zeros still
    // fall back to the established one-loop Standard-Model matching.
    if (sm_like_filter || (!bsm_only_generation && !full_target_generation)) {
        return MartyOrderPolicy::AUTO;
    }
    return MartyAdapter{}.get_marty_order_policy();
}

std::string order_policy_name(MartyOrderPolicy policy) {
    switch (policy) {
    case MartyOrderPolicy::TREE_LEVEL_ONLY:
        return "tree-level-only";
    case MartyOrderPolicy::ONE_LOOP_ONLY:
        return "one-loop-only";
    case MartyOrderPolicy::AUTO:
        return "auto";
    }
    return "auto";
}

std::string generation_mode_marker(const std::string& wilson,
                                   bool sm_like_filter,
                                   bool bsm_only_generation,
                                   bool full_target_generation,
                                   MartyOrderPolicy order_policy) {
    return "HYPERISO_MARTY_GENERATION_MODE: "
         + generation_mode(
             wilson,
             sm_like_filter,
             bsm_only_generation,
             full_target_generation
         )
         + "; order-policy=" + order_policy_name(order_policy);
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
        full_target_generation,
        effective_order_policy(sm_like_filter, bsm_split_generation, full_target_generation)
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
    std::shared_lock<std::shared_mutex> read_lock(marty_artifact_mutex);
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

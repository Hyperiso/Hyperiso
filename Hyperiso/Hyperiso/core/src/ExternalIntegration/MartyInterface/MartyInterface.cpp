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

constexpr const char* kMartyCacheAbi = "HYPERISO_MARTY_CACHE_ABI: pyhyperiso-1.0.4-v18";

constexpr const char* kMartyTreeRecipePrefix = "__HYPERISO_MARTY_TREE_RECIPE__|";
constexpr const char* kMartyTreeRecipeToken = "HYPERISO_MARTY_TREE_PROJECTION_TERMS";
constexpr const char* kMartyTreeRecipeInjectedBegin = "HYPERISO_MARTY_TREE_RECIPE_INJECTED_BEGIN";
constexpr const char* kMartyTreeRecipeInjectedEnd = "HYPERISO_MARTY_TREE_RECIPE_INJECTED_END";

struct MartyTreeProjectionTerm {
    std::string id;
    double weight {1.0};
    std::string left_current;
    std::string right_current;
    std::string layout;
    std::vector<int> fermion_order;
    std::vector<int> operator_order;
};

bool supports_tree_projection_recipe(const std::string& wilson) {
    if (wilson == "C9" || wilson == "C10"
        || wilson == "CP9" || wilson == "CP10") {
        return true;
    }
    // B -> s nu_i anti-nu_j uses the same four-fermion TreeLevel projection
    // machinery.  The CNU templates provide their own external-neutrino
    // insertions, while the runtime recipe controls F/O/current layout.
    return wilson.rfind("CNU_L_", 0) == 0 || wilson.rfind("CNU_R_", 0) == 0
        || wilson.rfind("CKNU_L_", 0) == 0 || wilson.rfind("CKNU_R_", 0) == 0;
}

std::vector<std::string> split_recipe_key(const std::string& key) {
    std::vector<std::string> fields;
    std::stringstream stream(key);
    std::string field;
    while (std::getline(stream, field, '|')) fields.push_back(field);
    return fields;
}

bool is_valid_marty_permutation(const std::vector<int>& order) {
    if (order.size() != 4) return false;
    auto sorted = order;
    std::sort(sorted.begin(), sorted.end());
    return sorted == std::vector<int>({0, 1, 2, 3});
}

std::string compact_order_marker(const std::vector<int>& order) {
    if (order.empty()) return "template-default";
    std::ostringstream out;
    for (std::size_t i = 0; i < order.size(); ++i) {
        if (i != 0) out << '-';
        out << order[i];
    }
    return out.str();
}

std::string cpp_dirac_coupling(const std::string& current) {
    if (current == "VL") return "mty::DiracCoupling::VL";
    if (current == "VR") return "mty::DiracCoupling::VR";
    if (current == "V")  return "mty::DiracCoupling::V";
    if (current == "A")  return "mty::DiracCoupling::A";
    throw std::runtime_error(
        "Unsupported MARTY tree projection current '" + current
        + "'. Recipe ABI v1 supports VL, VR, V and A."
    );
}

std::string cpp_vector_literal(const std::vector<int>& values) {
    std::ostringstream out;
    out << '{';
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i != 0) out << ", ";
        out << values[i];
    }
    out << '}';
    return out.str();
}

std::string escape_cpp_string(const std::string& value) {
    std::string out;
    out.reserve(value.size());
    for (const char c : value) {
        if (c == '\\' || c == '"') out.push_back('\\');
        out.push_back(c);
    }
    return out;
}

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
std::vector<MartyTreeProjectionTerm> effective_tree_projection_recipe(
    const std::string& wilson,
    bool sm_like_filter,
    bool bsm_only_generation,
    bool full_target_generation
);
std::string tree_projection_recipe_marker(const std::vector<MartyTreeProjectionTerm>& recipe);
void inject_tree_projection_recipe(const fs::path& generated_file,
                                   const std::vector<MartyTreeProjectionTerm>& recipe);
MartyOrderPolicy effective_order_policy(bool sm_like_filter,
                                          bool bsm_only_generation,
                                          bool full_target_generation);
std::vector<int> effective_fermion_order(const std::string& wilson,
                                         bool one_loop,
                                         bool sm_like_filter,
                                         bool bsm_only_generation,
                                         bool full_target_generation);
std::vector<int> effective_operator_order(const std::string& wilson,
                                          bool one_loop,
                                          bool sm_like_filter,
                                          bool bsm_only_generation,
                                          bool full_target_generation);
std::string generation_mode_marker(const std::string& wilson,
                                   bool sm_like_filter,
                                   bool bsm_only_generation,
                                   bool full_target_generation,
                                   MartyOrderPolicy order_policy,
                                   const std::vector<int>& tree_fermion_order,
                                   const std::vector<int>& one_loop_fermion_order,
                                   const std::vector<int>& tree_operator_order,
                                   const std::vector<int>& one_loop_operator_order);
void append_cache_metadata_if_missing(const fs::path& generated_file,
                                      const std::string& model_signature,
                                      const std::string& template_signature_value,
                                      const std::string& mode_marker);
bool template_needs_generic_tree_first(
    const std::string& wilson,
    const std::shared_ptr<FileNameManager>& files,
    bool bsm_only_generation,
    bool full_target_generation
);
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
    const bool tree_first_fallback = template_needs_generic_tree_first(
        wilson,
        files,
        bsm_split_generation,
        full_target_generation
    );
    const MartyOrderPolicy order_policy = effective_order_policy(
        sm_like_filter, bsm_split_generation, full_target_generation
    );
    const std::vector<int> tree_fermion_order = effective_fermion_order(
        wilson, false, sm_like_filter, bsm_split_generation, full_target_generation
    );
    const std::vector<int> one_loop_fermion_order = effective_fermion_order(
        wilson, true, sm_like_filter, bsm_split_generation, full_target_generation
    );
    const std::vector<int> tree_operator_order = effective_operator_order(
        wilson, false, sm_like_filter, bsm_split_generation, full_target_generation
    );
    const std::vector<int> one_loop_operator_order = effective_operator_order(
        wilson, true, sm_like_filter, bsm_split_generation, full_target_generation
    );
    invalidate_template_model_cache_if_needed(
        wilson, output_model, target_model, model_path, model_template_index,
        sm_like_filter, bsm_split_generation, full_target_generation
    );

    std::unique_ptr<ModelModifier> smModifier;
    smModifier = std::make_unique<GeneralModelModifier>(
        wilson, output_model, target_model, model_path, model_template_index,
        sm_like_filter, bsm_split_generation, full_target_generation,
        tree_first_fallback, order_policy,
        tree_fermion_order, one_loop_fermion_order,
        tree_operator_order, one_loop_operator_order
    );

    std::unique_ptr<TemplateManagerBase> templateManager = std::make_unique<NonNumericTemplateManager>(files->getTemplateDir());
    templateManager->setModelAndWilson(output_model, wilson);
    templateManager->setModelModifier(std::move(smModifier));

    CodeGenerator codeGenerator(std::move(templateManager));

    codeGenerator.generate(wilson, files->getGeneratedFileName());
    inject_tree_projection_recipe(
        files->getGeneratedFileName(),
        effective_tree_projection_recipe(
            wilson, sm_like_filter, bsm_split_generation, full_target_generation
        )
    );
    append_cache_metadata_if_missing(
        files->getGeneratedFileName(),
        GeneralModelModifier::modelSignature(target_model, model_path, model_template_index),
        template_signature(wilson, files),
        generation_mode_marker(
            wilson, sm_like_filter, bsm_split_generation, full_target_generation,
            order_policy,
            tree_fermion_order, one_loop_fermion_order,
            tree_operator_order, one_loop_operator_order
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
        effective_order_policy(sm_like_filter, bsm_split_generation, full_target_generation),
        effective_fermion_order(
            wilson, false, sm_like_filter, bsm_split_generation, full_target_generation
        ),
        effective_fermion_order(
            wilson, true, sm_like_filter, bsm_split_generation, full_target_generation
        ),
        effective_operator_order(
            wilson, false, sm_like_filter, bsm_split_generation, full_target_generation
        ),
        effective_operator_order(
            wilson, true, sm_like_filter, bsm_split_generation, full_target_generation
        )
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

    const auto execute_prepared_group = [&]() -> std::string {
        const auto prepared_it = prepared_groups_by_wilson.find(wilson);
        if (prepared_it == prepared_groups_by_wilson.end()) {
            return execute_isolated();
        }
        const PreparedGroup group = prepared_it->second;
        const fs::path output_root = FileNameManager::getInstance(wilson, output_model)->getOutputDir();
        const fs::path group_dir = output_root / "groups" /
            (sanitize_path_component(output_model) + "_" + sanitize_path_component(group.group));
        const fs::path cache_dir = group_dir / "numeric_cache";
        fs::create_directories(cache_dir);

        // Build a deterministic point signature from Q_match and every numeric
        // parameter consumed by every coefficient in the group. This is the
        // crucial guard against reusing a group result after switching LHA point.
        std::ostringstream signature;
        signature << kMartyCacheAbi << "|" << std::setprecision(17) << "Q=" << Q_match;
        for (const auto& member : group.members) {
            // Numeric point caches survive across processes.  Tie them to the
            // actual numeric executable as well as to the parameter values so a
            // rebuilt wrapper can never reuse a CSV produced by older code.
            const auto member_files = FileNameManager::getInstance(member, group.output_model);
            std::error_code executable_time_ec;
            const auto executable_time = fs::last_write_time(
                member_files->getNumExecutableFileName(), executable_time_ec
            );
            signature << "|NUMEXE=" << member << "@"
                      << (executable_time_ec
                              ? 0
                              : executable_time.time_since_epoch().count());
            auto params = snapshot_numeric_params(
                member, group.output_model, group.target_model,
                group.bsm_split_generation, group.full_target_generation
            );
            std::vector<std::pair<std::string, double>> ordered(params.begin(), params.end());
            std::sort(ordered.begin(), ordered.end(), [](const auto& lhs, const auto& rhs) {
                return lhs.first < rhs.first;
            });
            signature << "|" << member;
            for (const auto& [name, value] : ordered) {
                signature << ";" << name << "=" << std::setprecision(17) << value;
            }
        }
        const std::size_t point_hash = std::hash<std::string>{}(signature.str());
        const fs::path cache_file = cache_dir / ("point_" + std::to_string(point_hash) + ".csv");

        static std::mutex group_numeric_mutex;
        std::lock_guard<std::mutex> cache_lock(group_numeric_mutex);
        std::error_code ec;
        if (!fs::is_regular_file(cache_file, ec) || ec || fs::file_size(cache_file, ec) == 0 || ec) {
            fs::remove(cache_file, ec);
            for (const auto& member : group.members) {
                const auto member_files = FileNameManager::getInstance(member, group.output_model);
                const fs::path member_run_dir = make_invocation_directory(
                    member_files, member, group.output_model
                );
                const fs::path param_file = member_run_dir / "paramlist.csv";
                const fs::path output_file = member_run_dir / "wilson.csv";
                try {
                    write_parameter_snapshot(
                        param_file,
                        snapshot_numeric_params(
                            member, group.output_model, group.target_model,
                            group.bsm_split_generation, group.full_target_generation
                        )
                    );
                    MakeCompilerStrategy compiler(group.output_model, member);
                    compiler.set_Q_match(Q_match);
                    compiler.set_param_file(param_file);
                    compiler.set_output_file(output_file);
                    compiler.compile_run(
                        member_files->getLibDir(), member_files->getNumExecutableFileName()
                    );
                    publish_legacy_csv(output_file, cache_file);
                } catch (...) {
                    std::error_code cleanup_ec;
                    fs::remove_all(member_run_dir, cleanup_ec);
                    fs::remove(cache_file, cleanup_ec);
                    throw;
                }
                std::error_code cleanup_ec;
                fs::remove_all(member_run_dir, cleanup_ec);
            }
        }

        // Return an invocation-local copy so MartyWilson's existing cleanup
        // ownership remains correct. The persistent merged cache itself is never
        // handed to the coefficient object.
        const auto request_files = FileNameManager::getInstance(wilson, output_model);
        const fs::path run_dir = make_invocation_directory(request_files, wilson, output_model);
        const fs::path output_file = run_dir / "wilson.csv";
        fs::copy_file(cache_file, output_file, fs::copy_options::overwrite_existing, ec);
        if (ec) {
            fs::remove_all(run_dir, ec);
            throw std::runtime_error(
                "Cannot copy MARTY group numeric cache for " + wilson + ": " + ec.message()
            );
        }
        return output_file.string();
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
            return execute_prepared_group();
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
    return execute_prepared_group();
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

bool template_needs_generic_tree_first(
    const std::string& wilson,
    const std::shared_ptr<FileNameManager>& files,
    bool bsm_only_generation,
    bool full_target_generation
) {
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

    const bool has_tree_call = std::regex_search(source, tree_call);
    const bool has_loop_call = std::regex_search(source, loop_call);

    // Keep the historical tree probe for templates whose native leading call
    // is OneLoop.  For a pure BSM target, also wrap TreeLevel-only templates so
    // AUTO has one uniform meaning: evaluate LO first and fall back to NLO only
    // when the complete LO coefficient is structurally zero.  Templates that
    // already contain both orders (for example C10) retain their specialised
    // internal tree-first implementation.
    if (has_loop_call && !has_tree_call) {
        return true;
    }
    return bsm_only_generation
        && !full_target_generation
        && has_tree_call
        && !has_loop_call;
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
    if (sm_like_filter || !bsm_only_generation || full_target_generation) {
        return MartyOrderPolicy::AUTO;
    }
    return MartyAdapter{}.get_marty_order_policy();
}

std::vector<MartyTreeProjectionTerm> effective_tree_projection_recipe(
    const std::string& wilson,
    bool sm_like_filter,
    bool bsm_only_generation,
    bool full_target_generation
) {
    if (!supports_tree_projection_recipe(wilson)
        || sm_like_filter
        || (!bsm_only_generation && !full_target_generation)) {
        return {};
    }

    const MartyAdapter adapter;
    const auto fermion_orders = adapter.get_marty_tree_fermion_orders();
    const auto operator_orders = adapter.get_marty_tree_operator_orders();
    const std::string prefix = std::string(kMartyTreeRecipePrefix) + wilson + "|";

    std::vector<MartyTreeProjectionTerm> recipe;
    for (const auto& [key, fermion_order] : fermion_orders) {
        if (key.rfind(prefix, 0) != 0) continue;

        const auto fields = split_recipe_key(key);
        if (fields.size() != 7
            || fields[0] != "__HYPERISO_MARTY_TREE_RECIPE__"
            || fields[1] != wilson) {
            throw std::runtime_error("Malformed MARTY tree projection recipe key: " + key);
        }
        const auto operator_it = operator_orders.find(key);
        if (operator_it == operator_orders.end()) {
            throw std::runtime_error(
                "MARTY tree projection recipe term is missing its operator order: " + key
            );
        }
        if (!is_valid_marty_permutation(fermion_order)
            || !is_valid_marty_permutation(operator_it->second)) {
            throw std::runtime_error(
                "MARTY tree projection recipe F/O must be permutations of 0,1,2,3: " + key
            );
        }
        if (fields[6] != "quark_first" && fields[6] != "lepton_first") {
            throw std::runtime_error(
                "MARTY tree projection recipe layout must be quark_first or lepton_first: " + key
            );
        }
        (void)cpp_dirac_coupling(fields[4]);
        (void)cpp_dirac_coupling(fields[5]);

        MartyTreeProjectionTerm term;
        term.id = fields[2];
        term.weight = std::stod(fields[3]);
        if (!std::isfinite(term.weight)) {
            throw std::runtime_error(
                "MARTY tree projection recipe weight must be finite: " + key
            );
        }
        term.left_current = fields[4];
        term.right_current = fields[5];
        term.layout = fields[6];
        term.fermion_order = fermion_order;
        term.operator_order = operator_it->second;
        recipe.push_back(std::move(term));
    }

    for (const auto& [key, operator_order] : operator_orders) {
        if (key.rfind(prefix, 0) != 0) continue;
        if (fermion_orders.find(key) == fermion_orders.end()) {
            throw std::runtime_error(
                "MARTY tree projection recipe term is missing its fermion order: " + key
            );
        }
        (void)operator_order;
    }

    std::sort(recipe.begin(), recipe.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.id < rhs.id;
    });
    return recipe;
}

std::string tree_projection_recipe_marker(const std::vector<MartyTreeProjectionTerm>& recipe) {
    if (recipe.empty()) return "template-default-direct";
    std::ostringstream out;
    out << "recipe-v1[";
    for (std::size_t i = 0; i < recipe.size(); ++i) {
        const auto& term = recipe[i];
        if (i != 0) out << ';';
        out << term.id << ':' << std::setprecision(17) << term.weight
            << ':' << term.left_current << ',' << term.right_current
            << ':' << term.layout
            << ":F=" << compact_order_marker(term.fermion_order)
            << ":O=" << compact_order_marker(term.operator_order);
    }
    out << ']';
    return out.str();
}

void inject_tree_projection_recipe(
    const fs::path& generated_file,
    const std::vector<MartyTreeProjectionTerm>& recipe
) {
    std::ifstream input(generated_file);
    if (!input) {
        throw std::runtime_error(
            "Cannot inject MARTY tree projection recipe into: " + generated_file.string()
        );
    }
    std::string source(
        (std::istreambuf_iterator<char>(input)),
        std::istreambuf_iterator<char>()
    );
    const std::string token = kMartyTreeRecipeToken;
    const std::string begin_marker = std::string("// ") + kMartyTreeRecipeInjectedBegin;
    const std::string end_marker = std::string("// ") + kMartyTreeRecipeInjectedEnd;

    // Always keep a replaceable marker block in the generated source.  The old
    // v11/v12 implementation consumed the template token permanently; a second
    // build on the same cache then failed because there was no hook left to
    // inject into.  Marker-to-marker replacement is idempotent and also allows
    // a changed recipe to update an already generated source safely.
    std::ostringstream replacement;
    replacement << begin_marker << "\n";
    for (std::size_t i = 0; i < recipe.size(); ++i) {
        const auto& term = recipe[i];
        if (i != 0) replacement << ",\n";
        replacement
            << "        {\"" << escape_cpp_string(term.id) << "\", "
            << std::setprecision(17) << term.weight << ", "
            << cpp_vector_literal(term.fermion_order) << ", "
            << cpp_vector_literal(term.operator_order) << ", "
            << cpp_dirac_coupling(term.left_current) << ", "
            << cpp_dirac_coupling(term.right_current) << ", "
            << (term.layout == "lepton_first" ? "true" : "false") << "}";
    }
    if (!recipe.empty()) {
        replacement << "\n";
    }
    replacement << end_marker;
    const std::string replacement_text = replacement.str();

    std::size_t token_pos = source.find(token);
    if (token_pos != std::string::npos) {
        while ((token_pos = source.find(token, token_pos)) != std::string::npos) {
            source.replace(token_pos, token.size(), replacement_text);
            token_pos += replacement_text.size();
        }
    } else {
        const auto begin = source.find(begin_marker);
        const auto end = begin == std::string::npos
            ? std::string::npos
            : source.find(end_marker, begin + begin_marker.size());
        if (begin == std::string::npos || end == std::string::npos) {
            if (!recipe.empty()) {
                throw std::runtime_error(
                    "Configured MARTY tree projection recipe for a template without a recipe hook: "
                    + generated_file.string()
                );
            }
            return;
        }
        source.replace(
            begin,
            end + end_marker.size() - begin,
            replacement_text
        );
    }

    std::ofstream output(generated_file, std::ios::trunc);
    if (!output) {
        throw std::runtime_error(
            "Cannot write MARTY source after recipe injection: " + generated_file.string()
        );
    }
    output << source;
}

std::vector<int> effective_fermion_order(const std::string& wilson,
                                         bool one_loop,
                                         bool sm_like_filter,
                                         bool bsm_only_generation,
                                         bool full_target_generation) {
    if (sm_like_filter || (!bsm_only_generation && !full_target_generation)) return {};
    if (!one_loop && !effective_tree_projection_recipe(
            wilson, sm_like_filter, bsm_only_generation, full_target_generation
        ).empty()) {
        return {};
    }
    const MartyAdapter adapter;
    const auto orders = one_loop
        ? adapter.get_marty_one_loop_fermion_orders()
        : adapter.get_marty_tree_fermion_orders();
    const auto it = orders.find(wilson);
    return it == orders.end() ? std::vector<int>{} : it->second;
}

std::vector<int> effective_operator_order(const std::string& wilson,
                                          bool one_loop,
                                          bool sm_like_filter,
                                          bool bsm_only_generation,
                                          bool full_target_generation) {
    if (sm_like_filter || (!bsm_only_generation && !full_target_generation)) return {};
    if (!one_loop && !effective_tree_projection_recipe(
            wilson, sm_like_filter, bsm_only_generation, full_target_generation
        ).empty()) {
        return {};
    }
    const MartyAdapter adapter;
    const auto orders = one_loop
        ? adapter.get_marty_one_loop_operator_orders()
        : adapter.get_marty_tree_operator_orders();
    const auto it = orders.find(wilson);
    return it == orders.end() ? std::vector<int>{} : it->second;
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

std::string fermion_order_marker(const std::vector<int>& fermion_order) {
    if (fermion_order.empty()) {
        return "template-default";
    }
    std::ostringstream order;
    for (std::size_t i = 0; i < fermion_order.size(); ++i) {
        if (i != 0) {
            order << '-';
        }
        order << fermion_order[i];
    }
    return order.str();
}

std::string generation_mode_marker(const std::string& wilson,
                                   bool sm_like_filter,
                                   bool bsm_only_generation,
                                   bool full_target_generation,
                                   MartyOrderPolicy order_policy,
                                   const std::vector<int>& tree_fermion_order,
                                   const std::vector<int>& one_loop_fermion_order,
                                   const std::vector<int>& tree_operator_order,
                                   const std::vector<int>& one_loop_operator_order) {
    return "HYPERISO_MARTY_GENERATION_MODE: "
         + generation_mode(
             wilson,
             sm_like_filter,
             bsm_only_generation,
             full_target_generation
         )
         + "; order-policy=" + order_policy_name(order_policy)
         + "; tree-fermion-order=" + fermion_order_marker(tree_fermion_order)
         + "; one-loop-fermion-order=" + fermion_order_marker(one_loop_fermion_order)
         + "; tree-operator-order=" + fermion_order_marker(tree_operator_order)
         + "; one-loop-operator-order=" + fermion_order_marker(one_loop_operator_order)
         + (supports_tree_projection_recipe(wilson)
                ? "; tree-projection=" + tree_projection_recipe_marker(
                    effective_tree_projection_recipe(
                        wilson, sm_like_filter, bsm_only_generation, full_target_generation
                    )
                  )
                : "");
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
        effective_order_policy(sm_like_filter, bsm_split_generation, full_target_generation),
        effective_fermion_order(
            wilson, false, sm_like_filter, bsm_split_generation, full_target_generation
        ),
        effective_fermion_order(
            wilson, true, sm_like_filter, bsm_split_generation, full_target_generation
        ),
        effective_operator_order(
            wilson, false, sm_like_filter, bsm_split_generation, full_target_generation
        ),
        effective_operator_order(
            wilson, true, sm_like_filter, bsm_split_generation, full_target_generation
        )
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


namespace {
bool hyperiso_nonempty_file(const fs::path& path) {
    std::error_code ec;
    return fs::is_regular_file(path, ec) && fs::file_size(path, ec) > 0;
}

bool hyperiso_file_is_at_least_as_new_as(const fs::path& candidate, const fs::path& dependency) {
    std::error_code ec;
    if (!hyperiso_nonempty_file(candidate) || !fs::exists(dependency, ec)) return false;
    const auto ctime = fs::last_write_time(candidate, ec);
    if (ec) return false;
    const auto dtime = fs::last_write_time(dependency, ec);
    if (ec) return false;
    return ctime >= dtime;
}

std::string hyperiso_read_text_if_exists(const fs::path& path) {
    std::ifstream input(path);
    if (!input) return {};
    return std::string((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
}

void hyperiso_write_text_if_changed(const fs::path& path, const std::string& content) {
    if (hyperiso_read_text_if_exists(path) == content) return;
    std::ofstream output(path, std::ios::trunc);
    if (!output) throw std::runtime_error("Cannot write MARTY group source: " + path.string());
    output << content;
}

std::string hyperiso_cpp_quote(const std::string& value) {
    std::string out;
    out.reserve(value.size() + 8);
    for (char c : value) {
        if (c == '\\' || c == '"') out.push_back('\\');
        out.push_back(c);
    }
    return out;
}


std::string hyperiso_plugin_source(const fs::path& source_path,
                                   const std::string& wilson,
                                   const std::string& model_instantiation,
                                   bool expected_nonzero,
                                   bool scan_allowed,
                                   bool explicit_recipe) {
    std::ifstream input(source_path);
    if (!input) throw std::runtime_error("Cannot read generated MARTY source: " + source_path.string());
    std::string source((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());

    // Track the primary coefficient expression without knowing the local variable
    // name used by the individual template or specialised generated main().
    const std::regex add_re(
        "wilsonLib\\.addFunction\\(\\\"" + wilson + "\\\"\\s*,\\s*([^;]+)\\);"
    );
    std::smatch add_match;
    if (!std::regex_search(source, add_match, add_re)) {
        throw std::runtime_error("MARTY group batching cannot identify wilsonLib.addFunction for " + wilson);
    }
    const std::string expression = add_match[1].str();
    const std::string add_line = add_match[0].str();
    const std::string checked_add =
        "hyperiso_marty_group_primary_nonzero = (DeepRefreshed(" + expression + ") != CSL_0);\n    " + add_line;
    source.replace(add_match.position(0), add_match.length(0), checked_add);

    const std::string declaration =
        "\nstatic bool hyperiso_marty_group_primary_nonzero = true;\n";
    const auto using_pos = source.find("using namespace");
    if (using_pos == std::string::npos) throw std::runtime_error("Cannot inject MARTY group state for " + wilson);
    source.insert(using_pos, declaration);

    // GeneralModelModifier may replace the tiny template main() with a much
    // larger Tree-first or reg_prop main(). Transform that generated main as a
    // whole rather than relying on a particular calculate_* function name.
    const std::regex main_signature(R"(int\s+main\s*\(\s*\)\s*\{)");
    std::smatch main_match;
    if (!std::regex_search(source, main_match, main_signature)) {
        throw std::runtime_error("MARTY group batching cannot find generated main() for " + wilson);
    }
    const std::size_t main_begin = static_cast<std::size_t>(main_match.position(0));
    const std::size_t open_brace = source.find('{', main_begin);
    std::size_t cursor = open_brace + 1;
    int depth = 1;
    bool in_string = false;
    bool in_char = false;
    bool escape = false;
    for (; cursor < source.size() && depth > 0; ++cursor) {
        const char c = source[cursor];
        if (escape) { escape = false; continue; }
        if ((in_string || in_char) && c == '\\') { escape = true; continue; }
        if (!in_char && c == '"') { in_string = !in_string; continue; }
        if (!in_string && c == '\'') { in_char = !in_char; continue; }
        if (in_string || in_char) continue;
        if (c == '{') ++depth;
        else if (c == '}') --depth;
    }
    if (depth != 0) {
        throw std::runtime_error("MARTY group batching found an unterminated main() for " + wilson);
    }
    const std::size_t main_end = cursor; // one past the closing brace
    std::string body = source.substr(open_brace + 1, main_end - open_brace - 2);

    // Replace exactly one primary target-model construction. For C9/AUTO this
    // is tree_model; its deliberately isolated loop_model fallback remains a
    // fresh target-model instance and therefore preserves the existing MARTY
    // global-state safety rule.
    const std::vector<std::string> shared_names = {"tree_model", "model", "sm"};
    bool shared_model_installed = false;
    for (const auto& name : shared_names) {
        const std::string declaration_line = model_instantiation + " " + name + ";";
        const auto pos = body.find(declaration_line);
        if (pos != std::string::npos) {
            body.replace(pos, declaration_line.size(), "mty::Model& " + name + " = *hyperiso_marty_shared_model;");
            shared_model_installed = true;
            break;
        }
    }
    if (!shared_model_installed) {
        throw std::runtime_error(
            "MARTY group batching cannot identify the primary target-model construction for " + wilson
        );
    }

    // The untouched simple-template main returns calculate_*(...) directly.
    // Generated Tree-first mains instead have a final `return 0;`. Support both.
    static const std::regex direct_return(
        R"(return\s+([A-Za-z_][A-Za-z0-9_]*\s*\([^;]+\))\s*;)"
    );
    std::smatch direct_match;
    if (std::regex_search(body, direct_match, direct_return)) {
        const std::string call = direct_match[1].str();
        const std::string replacement =
            "const int hyperiso_marty_group_status = " + call + ";\n"
            "    if (hyperiso_marty_group_status != 0) return hyperiso_marty_group_status;\n"
            "    return hyperiso_marty_group_primary_nonzero ? 0 : 42;";
        body.replace(direct_match.position(0), direct_match.length(0), replacement);
    } else {
        const std::string return_zero = "return 0;";
        const auto return_pos = body.rfind(return_zero);
        if (return_pos == std::string::npos) {
            throw std::runtime_error("MARTY group batching cannot identify main() return for " + wilson);
        }
        body.replace(return_pos, return_zero.size(),
                     "return hyperiso_marty_group_primary_nonzero ? 0 : 42;");
    }

    std::ostringstream replacement;
    replacement << "extern \"C\" int hyperiso_marty_group_entry(mty::Model* hyperiso_marty_shared_model) {\n"
                << "    if (hyperiso_marty_shared_model == nullptr) return 90;\n"
                << "    mty::Model::current = hyperiso_marty_shared_model;\n"
                << "    hyperiso_marty_group_primary_nonzero = true;\n"
                << body << "\n}\n"
                << "extern \"C\" int hyperiso_marty_group_expected() { return " << (expected_nonzero ? 1 : 0) << "; }\n"
                << "extern \"C\" int hyperiso_marty_group_scan_allowed() { return " << (scan_allowed ? 1 : 0) << "; }\n"
                << "extern \"C\" int hyperiso_marty_group_explicit_recipe() { return " << (explicit_recipe ? 1 : 0) << "; }\n";
    source.replace(main_begin, main_end - main_begin, replacement.str());
    return source;
}
} // namespace

bool MartyInterface::prepare_group(const std::string& group,
                                   const std::vector<std::string>& input_members,
                                   const std::string& output_model,
                                   const std::string& target_model,
                                   const std::string& model_path,
                                   bool sm_like_filter,
                                   bool bsm_split_generation,
                                   bool full_target_generation) {
    const MartyAdapter adapter;
    // A new preparation attempt supersedes any previous in-process group state.
    // This is important when users switch batching/expected-nonzero settings and
    // rebuild without restarting Python.
    for (const auto& member : input_members) prepared_groups_by_wilson.erase(member);
    if (!adapter.get_marty_group_batching() || input_members.empty()) return false;
    if (!MartyRuntimeConfig::require_available("MartyInterface::prepare_group").valid) return false;

    std::vector<std::string> members = input_members;
    // C9/CP9/CP10 have specialised reg_prop/global-state handling. Run them
    // after ordinary coefficients, with C9 last, so a fallback cannot affect a
    // later plugin using the shared model.
    std::stable_sort(members.begin(), members.end(), [](const std::string& a, const std::string& b) {
        auto rank = [](const std::string& x) {
            if (x == "C9") return 2;
            if (x == "CP9" || x == "CP10") return 1;
            return 0;
        };
        return rank(a) < rank(b);
    });

    const auto expected_vector = adapter.get_marty_expected_nonzero_coefficients();
    const std::unordered_set<std::string> expected(expected_vector.begin(), expected_vector.end());

    // Refuse a partially batchable group. Silent hybrid analytical state is much
    // harder to reason about than a clean fallback to the historical path.
    for (const auto& wilson : members) {
        const auto files = FileNameManager::getInstance(wilson, output_model);
        if (!fs::is_regular_file(fs::path(files->getTemplateDir()) / (wilson + ".cpp"))) {
            LOG_WARN("MartyGroupBatch", "Group ", group, " contains no MARTY template for ", wilson,
                     "; falling back to coefficient-by-coefficient generation.");
            return false;
        }
    }

    const auto model_template_index = resolve_model_template_index(target_model);
    const std::string model_instantiation = GeneralModelModifier::resolveModelInstantiation(
        target_model, model_path, model_template_index
    );
    const fs::path output_root = FileNameManager::getInstance(members.front(), output_model)->getOutputDir();
    const fs::path group_dir = output_root / "groups" /
        (sanitize_path_component(output_model) + "_" + sanitize_path_component(group));
    fs::create_directories(group_dir);

    // Persistent analytical group signature.  A successful group build is a
    // model/template/order/projection artifact, not a parameter-point artifact:
    // once it exists it can be reused by later Python processes.  Numeric point
    // caching remains separate and is keyed by Q_match + the actual parameters.
    std::ostringstream group_signature_stream;
    group_signature_stream << kMartyCacheAbi << "\n"
                           << "HYPERISO_MARTY_GROUP_CACHE_SCHEMA: v1\n"
                           << "group=" << group << "\n"
                           << "output_model=" << output_model << "\n"
                           << "target_model=" << target_model << "\n"
                           << "model_instantiation=" << model_instantiation << "\n"
                           << "model_signature="
                           << GeneralModelModifier::modelSignature(
                                  target_model, model_path, model_template_index
                              )
                           << "\n"
                           << "sm_like_filter=" << (sm_like_filter ? 1 : 0) << "\n"
                           << "bsm_split_generation=" << (bsm_split_generation ? 1 : 0) << "\n"
                           << "full_target_generation=" << (full_target_generation ? 1 : 0) << "\n";
    {
        std::vector<std::string> expected_sorted(expected_vector.begin(), expected_vector.end());
        std::sort(expected_sorted.begin(), expected_sorted.end());
        group_signature_stream << "expected_nonzero=";
        for (const auto& name : expected_sorted) group_signature_stream << name << ",";
        group_signature_stream << "\n";
    }
    for (const auto& wilson : members) {
        const auto files = FileNameManager::getInstance(wilson, output_model);
        group_signature_stream
            << "member=" << wilson << "\n"
            << template_signature(wilson, files) << "\n"
            << generation_mode_marker(
                   wilson, sm_like_filter, bsm_split_generation, full_target_generation,
                   effective_order_policy(
                       sm_like_filter, bsm_split_generation, full_target_generation
                   ),
                   effective_fermion_order(
                       wilson, false, sm_like_filter, bsm_split_generation,
                       full_target_generation
                   ),
                   effective_fermion_order(
                       wilson, true, sm_like_filter, bsm_split_generation,
                       full_target_generation
                   ),
                   effective_operator_order(
                       wilson, false, sm_like_filter, bsm_split_generation,
                       full_target_generation
                   ),
                   effective_operator_order(
                       wilson, true, sm_like_filter, bsm_split_generation,
                       full_target_generation
                   )
               )
            << "\n";
    }
    const std::string group_signature = group_signature_stream.str();
    const fs::path analytical_ready_file = group_dir / "analytical.ready";

    const auto non_empty_file = [](const fs::path& path) {
        std::error_code ec;
        return fs::is_regular_file(path, ec) && !ec
            && fs::file_size(path, ec) > 0 && !ec;
    };
    const auto persistent_group_artifacts_ready = [&]() {
        if (hyperiso_read_text_if_exists(analytical_ready_file) != group_signature) return false;
        for (const auto& wilson : members) {
            const auto files = FileNameManager::getInstance(wilson, output_model);
            const fs::path lib_dir = files->getLibDir();
            if (!non_empty_file(files->getGeneratedFileName())
                || !non_empty_file(files->getExecutableFileName())
                || !non_empty_file(files->getNumGeneratedFileName())
                || !non_empty_file(files->getNumExecutableFileName())
                || !fs::is_regular_file(lib_dir / "Makefile")) {
                return false;
            }
        }
        return true;
    };

    const auto register_prepared_group = [&]() {
        PreparedGroup prepared{group, members, output_model, target_model, model_path,
                               sm_like_filter, bsm_split_generation, full_target_generation};
        for (const auto& wilson : members) prepared_groups_by_wilson[wilson] = prepared;
    };

    // On a cache hit there is no analytical work at all: no plugin generation,
    // no model construction and no MARTY process/matching.  In a fresh Python
    // process we still regenerate the lightweight numeric wrapper metadata once
    // to repopulate the in-memory dependency map; its executable is reused when
    // unchanged.  In the same process even that step is skipped.
    if (persistent_group_artifacts_ready()) {
        std::cout << "[MARTY group " << group
                  << "] analytical cache HIT: reusing validated MARTY libraries; "
                     "skipping stages 1-3."
                  << std::endl;
        for (const auto& wilson : members) {
            if (!dependencies.contains(wilson)) {
                generate_numlib(wilson, output_model, target_model,
                                bsm_split_generation, full_target_generation);
                compile_numlib(wilson, output_model);
            }
        }
        register_prepared_group();
        return true;
    }

    std::cout << "[MARTY group " << group
              << "] analytical cache MISS: rebuilding group artifacts." << std::endl;
    // Only an analytical cache miss invalidates point-wise numerical caches.
    {
        std::error_code ec;
        fs::remove_all(group_dir / "numeric_cache", ec);
        fs::remove(analytical_ready_file, ec);
    }

    struct Plugin {
        std::string wilson;
        fs::path path;
        bool four_fermion {false};
        bool expected_nonzero {false};
        bool scan_allowed {false};
        bool explicit_recipe {false};
    };
    std::vector<Plugin> plugins;
    plugins.reserve(members.size());

    const auto group_prepare_start = std::chrono::steady_clock::now();
    std::cout << "[MARTY group " << group << "] analytical batch:\n"
              << members.size() << " coefficients, shared base-model instance." << std::endl;
    std::cout << "[MARTY group " << group << "] stage 1/3: preparing analytical plugins"
              << std::endl;

    try {
        std::size_t plugin_index = 0;
        for (const auto& wilson : members) {
            ++plugin_index;
            const auto plugin_prepare_start = std::chrono::steady_clock::now();
            std::cout << "[MARTY group " << group << "] [plugin " << plugin_index << "/"
                      << members.size() << "] " << wilson << ": generating source ..."
                      << std::endl;
            generate(wilson, output_model, target_model, model_path, sm_like_filter,
                     bsm_split_generation, full_target_generation);
            const auto files = FileNameManager::getInstance(wilson, output_model);
            std::ifstream generated(files->getGeneratedFileName());
            std::string generated_source((std::istreambuf_iterator<char>(generated)), std::istreambuf_iterator<char>());

            // Determine scan eligibility from the coefficient template itself, not
            // from the generated source.  GeneralModelModifier injects the
            // hyperiso_marty_dimension6_operator helper into generated sources even
            // for dipoles/custom projectors (C7/C8, C5/C6), which previously made
            // those coefficients look F/O-scanable although their actual projector
            // does not use the generic dimension-6 operator-order machinery.
            const fs::path coefficient_template =
                fs::path(files->getTemplateDir()) / (wilson + ".cpp");
            std::ifstream template_input(coefficient_template);
            std::string template_source((std::istreambuf_iterator<char>(template_input)),
                                        std::istreambuf_iterator<char>());
            const bool four_fermion =
                template_source.find("dimension6Operator(") != std::string::npos
                || template_source.find("hyperiso_marty_dimension6_operator(") != std::string::npos;
            const bool explicit_recipe = generated_source.find("recipe-v1[") != std::string::npos;
            const bool scan_allowed = four_fermion && !explicit_recipe;

            const fs::path plugin_cpp = group_dir / (sanitize_path_component(wilson) + "_plugin.cpp");
            const fs::path plugin_so = group_dir / ("lib" + sanitize_path_component(wilson) + ".so");
            const std::string plugin_source = hyperiso_plugin_source(
                files->getGeneratedFileName(), wilson, model_instantiation,
                expected.contains(wilson), scan_allowed, explicit_recipe
            );
            hyperiso_write_text_if_changed(plugin_cpp, plugin_source);
            std::cout << "[MARTY group " << group << "] [plugin " << plugin_index << "/"
                      << members.size() << "] " << wilson
                      << ": source ready | projector=" << (four_fermion ? "dimension6" : "template-specific")
                      << " | expected=" << (expected.contains(wilson) ? "nonzero" : "optional")
                      << " | recipe=" << (explicit_recipe ? "explicit" : "none")
                      << " | scan=" << (scan_allowed ? "available" : "disabled")
                      << std::endl;

            const auto active_marty = MartyRuntimeConfig::require_available("MartyInterface::prepare_group plugin cache");
            const bool plugin_cached = hyperiso_file_is_at_least_as_new_as(plugin_so, plugin_cpp)
                && (!active_marty.valid
                    || hyperiso_file_is_at_least_as_new_as(plugin_so, active_marty.marty_library));
            if (plugin_cached) {
                std::cout << "[MARTY group " << group << "] [plugin " << plugin_index << "/"
                          << members.size() << "] " << wilson << ": cached shared object reused"
                          << std::endl;
            } else {
                std::cout << "[MARTY group " << group << "] [plugin " << plugin_index << "/"
                          << members.size() << "] " << wilson << ": compiling plugin ..."
                          << std::endl;
                GppCompilerStrategy compiler(output_model, wilson);
                compiler.compile_shared(plugin_cpp.string(), plugin_so.string());
            }
            const auto plugin_prepare_stop = std::chrono::steady_clock::now();
            std::cout << "[MARTY group " << group << "] [plugin " << plugin_index << "/"
                      << members.size() << "] " << wilson << ": plugin ready in "
                      << std::chrono::duration<double>(plugin_prepare_stop - plugin_prepare_start).count()
                      << " s" << std::endl;
            plugins.push_back({wilson, plugin_so, four_fermion, expected.contains(wilson),
                               scan_allowed, explicit_recipe});
        }

        std::cout << "[MARTY group " << group
                  << "] stage 1/3 complete: all plugins ready; launching shared-model driver."
                  << std::endl;

        const fs::path driver_cpp = group_dir / "group_driver.cpp";
        const fs::path driver_bin = group_dir / "group_driver";
        std::ofstream driver(driver_cpp, std::ios::trunc);
        driver << "#include <marty.h>\n#include <dlfcn.h>\n#include <array>\n#include <algorithm>\n"
               << "#include <chrono>\n#include <cstdlib>\n#include <iomanip>\n#include <iostream>\n#include <string>\n#include <vector>\n"
               << "#include \"" << hyperiso_cpp_quote(model_path) << "\"\n"
               << "using namespace mty; using namespace csl;\n"
               << "static std::string ord(const std::array<int,4>& a){ return std::to_string(a[0])+\",\"+std::to_string(a[1])+\",\"+std::to_string(a[2])+\",\"+std::to_string(a[3]); }\n"
               << "static double elapsed(std::chrono::steady_clock::time_point t){ return std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count(); }\n"
               << "int main(){ auto group_start=std::chrono::steady_clock::now(); std::cout<<\"[MARTY group " << hyperiso_cpp_quote(group) << "] stage 2/3: shared model construction START | model=" << hyperiso_cpp_quote(model_instantiation) << "\"<<std::endl; auto model_start=std::chrono::steady_clock::now(); mty::sm_input::undefineNumericalValues(); " << model_instantiation << " model; mty::Model::current=&model; std::cout<<\"[MARTY group " << hyperiso_cpp_quote(group) << "] stage 2/3: shared model READY | elapsed=\"<<elapsed(model_start)<<\" s\"<<std::endl;\n"
               << "std::cout<<\"[MARTY group " << hyperiso_cpp_quote(group) << "] stage 3/3: coefficient matching START | coefficients=" << members.size() << "\"<<std::endl;\n"
               << "std::array<int,4> base{0,1,2,3}; std::vector<std::array<int,4>> perms; do{perms.push_back(base);}while(std::next_permutation(base.begin(),base.end()));\n"
               << "int passed=0, expn=0, zero_allowed=0, failed=0;\n";
        std::size_t coefficient_index = 0;
        for (const auto& plugin : plugins) {
            ++coefficient_index;
            driver << "{ auto coefficient_start=std::chrono::steady_clock::now(); "
                   << "std::cout<<\"[MARTY group " << hyperiso_cpp_quote(group) << "] [coefficient "
                   << coefficient_index << "/" << members.size() << "] " << plugin.wilson
                   << " START | projector=" << (plugin.four_fermion ? "dimension6" : "template-specific")
                   << " | expected=" << (plugin.expected_nonzero ? "nonzero" : "optional")
                   << " | recipe=" << (plugin.explicit_recipe ? "explicit" : "none")
                   << " | scan=" << (plugin.scan_allowed ? "available" : "disabled")
                   << "\"<<std::endl; "
                   << "const char* path=\"" << hyperiso_cpp_quote(plugin.path.string()) << "\"; void* h=dlopen(path,RTLD_NOW|RTLD_LOCAL); if(!h){std::cerr<<dlerror()<<std::endl; return 91;}\n"
                   << "auto entry=reinterpret_cast<int(*)(mty::Model*)>(dlsym(h,\"hyperiso_marty_group_entry\"));"
                   << "auto expected=reinterpret_cast<int(*)()>(dlsym(h,\"hyperiso_marty_group_expected\"));"
                   << "auto scan=reinterpret_cast<int(*)()>(dlsym(h,\"hyperiso_marty_group_scan_allowed\"));"
                   << "auto recipe=reinterpret_cast<int(*)()>(dlsym(h,\"hyperiso_marty_group_explicit_recipe\"));"
                   << "if(!entry||!expected||!scan||!recipe) return 92; unsetenv(\"HYPERISO_MARTY_RUNTIME_TREE_F\"); unsetenv(\"HYPERISO_MARTY_RUNTIME_TREE_O\"); unsetenv(\"HYPERISO_MARTY_SCAN_QUIET\"); const bool is_expected=(expected()!=0); int st=entry(&model); if(is_expected)++expn; if(st!=0 && st!=42){std::cerr<<\"[MARTY group] " << plugin.wilson << " ERROR status=\"<<st<<std::endl; return st;} bool found=(st==0); bool zero=(st==42);\n"
                   << "if(zero && is_expected && scan() && !recipe()){ std::cout<<\"[MARTY group scan] " << plugin.wilson << ": configured/default projection is zero; scanning 24 F orders\"<<std::endl; setenv(\"HYPERISO_MARTY_SCAN_QUIET\",\"1\",1); "
                   << "for(std::size_t fi=0;fi<perms.size();++fi){const auto& f=perms[fi]; setenv(\"HYPERISO_MARTY_RUNTIME_TREE_F\",ord(f).c_str(),1); unsetenv(\"HYPERISO_MARTY_RUNTIME_TREE_O\"); st=entry(&model); if(st!=0 && st!=42) return st; if(st==0){found=true; std::cout<<\"[MARTY group scan] " << plugin.wilson << ": FOUND F=\"<<ord(f)<<\" O=configured/template\"<<std::endl; break;}}"
                   << "if(!found){std::cout<<\"[MARTY group scan] " << plugin.wilson << ": F scan exhausted; scanning 24 O orders\"<<std::endl; unsetenv(\"HYPERISO_MARTY_RUNTIME_TREE_F\"); for(std::size_t oi=0;oi<perms.size();++oi){const auto& o=perms[oi]; setenv(\"HYPERISO_MARTY_RUNTIME_TREE_O\",ord(o).c_str(),1); st=entry(&model); if(st!=0 && st!=42) return st; if(st==0){found=true; std::cout<<\"[MARTY group scan] " << plugin.wilson << ": FOUND F=configured/template O=\"<<ord(o)<<std::endl; break;}}}"
                   << "if(!found){std::cout<<\"[MARTY group scan] " << plugin.wilson << ": O scan exhausted; scanning F x O (576 maximum)\"<<std::endl; for(std::size_t fi=0;fi<perms.size();++fi){const auto& f=perms[fi]; setenv(\"HYPERISO_MARTY_RUNTIME_TREE_F\",ord(f).c_str(),1); if(fi==0 || (fi+1)%4==0) std::cout<<\"[MARTY group scan] " << plugin.wilson << ": F x O row \"<<(fi+1)<<\"/24 | F=\"<<ord(f)<<std::endl; for(const auto& o:perms){setenv(\"HYPERISO_MARTY_RUNTIME_TREE_O\",ord(o).c_str(),1); st=entry(&model); if(st!=0 && st!=42) return st; if(st==0){found=true; std::cout<<\"[MARTY group scan] " << plugin.wilson << ": FOUND F=\"<<ord(f)<<\" O=\"<<ord(o)<<std::endl; break;}} if(found)break;}} unsetenv(\"HYPERISO_MARTY_SCAN_QUIET\"); }\n"
                   << "unsetenv(\"HYPERISO_MARTY_RUNTIME_TREE_F\"); unsetenv(\"HYPERISO_MARTY_RUNTIME_TREE_O\"); unsetenv(\"HYPERISO_MARTY_SCAN_QUIET\");"
                   << "if(found){std::cout<<\"[MARTY group] " << plugin.wilson << "\"<<(is_expected?\" expected=nonzero\":\"\")<<\" PASS | elapsed=\"<<elapsed(coefficient_start)<<\" s\"<<std::endl; ++passed;} else if(!is_expected){std::cout<<\"[MARTY group] " << plugin.wilson << " ZERO (allowed) | elapsed=\"<<elapsed(coefficient_start)<<\" s\"<<std::endl; ++passed; ++zero_allowed;} else {std::cout<<\"[MARTY group] " << plugin.wilson << " expected=nonzero FAILED | elapsed=\"<<elapsed(coefficient_start)<<\" s\"<<std::endl; ++failed; return 42;} /* Keep plugin loaded until process exit: the shared MARTY model may retain objects/typeinfo from this TU. */ }\n";
        }
        driver << "std::cout<<\"[MARTY group " << hyperiso_cpp_quote(group) << "] summary: passed=\"<<passed<<\"/" << members.size() << ", expected-nonzero=\"<<expn<<\", zero-allowed=\"<<zero_allowed<<\", failed=\"<<failed<<\", analytical-elapsed=\"<<elapsed(group_start)<<\" s\"<<std::endl; return failed?1:0;}\n";
        driver.close();

        GppCompilerStrategy group_compiler(output_model, members.front());
        group_compiler.compile_group_driver(driver_cpp.string(), driver_bin.string());
        const std::string driver_command =
            "cd " + MartyRuntimeConfig::shell_quote(output_root)
            + " && " + MartyRuntimeConfig::shell_quote(driver_bin);
        if (!executeCommandStreaming(driver_command)) {
            throw std::runtime_error(
                "MARTY group analytical batch failed for " + group
                + "; numerical generation was not started."
            );
        }
        const auto group_prepare_stop = std::chrono::steady_clock::now();
        std::cout << "[MARTY group " << group
                  << "] analytical generation complete | host+driver elapsed="
                  << std::chrono::duration<double>(group_prepare_stop - group_prepare_start).count()
                  << " s; preparing numerical wrappers." << std::endl;

        // Each plugin has now emitted its numerical library while sharing the
        // same model. Build the lightweight numeric wrappers and dependencies.
        for (const auto& wilson : members) {
            generate_numlib(wilson, output_model, target_model,
                            bsm_split_generation, full_target_generation);
            compile_numlib(wilson, output_model);
            const auto files = FileNameManager::getInstance(wilson, output_model);
            std::ofstream marker(files->getExecutableFileName(), std::ios::trunc);
            marker << "HYPERISO_MARTY_GROUP_BATCHED\n" << group << "\n";
        }
        // Commit the persistent group cache only after both analytical libraries
        // and every numerical wrapper have completed successfully.
        hyperiso_write_text_if_changed(analytical_ready_file, group_signature);
    } catch (const std::exception& error) {
        for (const auto& wilson : members) prepared_groups_by_wilson.erase(wilson);
        throw;
    }

    register_prepared_group();
    return true;
}

bool MartyInterface::is_group_prepared(const std::string& wilson,
                                       const std::string& output_model,
                                       const std::string& target_model,
                                       const std::string& model_path) const {
    const auto it = prepared_groups_by_wilson.find(wilson);
    if (it == prepared_groups_by_wilson.end()) return false;
    const auto& group = it->second;
    return group.output_model == output_model
        && group.target_model == target_model
        && normalized_path(group.model_path) == normalized_path(model_path);
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

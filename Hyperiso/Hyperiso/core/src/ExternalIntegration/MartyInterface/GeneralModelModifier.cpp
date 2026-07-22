#include "GeneralModelModifier.h"
#include "MartyRuntimeConfig.h"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <utility>

namespace fs = std::filesystem;

namespace {

std::string stable_file_fingerprint(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("Cannot fingerprint MARTY model file: " + path.string());
    }

    // FNV-1a is intentionally simple and deterministic.  This is a cache key,
    // not a cryptographic integrity check.
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

} // namespace

GeneralModelModifier::GeneralModelModifier(std::string wilson,
                                             std::string model,
                                             std::string model_path,
                                             std::optional<int> model_template_index)
    : GeneralModelModifier(
        std::move(wilson),
        model,
        model,
        std::move(model_path),
        model_template_index,
        false,
        false,
        false,
        false,
        MartyOrderPolicy::AUTO,
        {}
      ) {}

GeneralModelModifier::GeneralModelModifier(std::string wilson,
                                             std::string output_model,
                                             std::string target_model,
                                             std::string model_path,
                                             std::optional<int> model_template_index,
                                             bool disable_non_sm_particles,
                                             bool bsm_split_generation,
                                             bool full_target_generation,
                                             bool tree_first_fallback,
                                             MartyOrderPolicy order_policy,
                                             std::vector<int> tree_fermion_order) {
        this->wilson = std::move(wilson);
        this->output_model = std::move(output_model);
        this->target_model = std::move(target_model);
        this->model_path = std::move(model_path);
        this->model_template_index = model_template_index;
        this->disable_non_sm_particles = disable_non_sm_particles;
        this->bsm_split_generation = bsm_split_generation;
        this->full_target_generation = full_target_generation;
        this->tree_first_fallback = tree_first_fallback;
        this->order_policy = order_policy;
        this->tree_fermion_order = std::move(tree_fermion_order);

        ModelFileChecker checker(this->model_path);
        this->model_class = checker.resolveModelClass(this->target_model);
        this->model_instantiation = resolveModelInstantiation(
            this->target_model,
            this->model_path,
            this->model_template_index
        );

        const auto marty = MartyRuntimeConfig::require_available("GeneralModelModifier");
        if (marty.valid) {
            this->marty_path = marty.marty_header.string();
        }
    }

std::string GeneralModelModifier::resolveModelInstantiation(const std::string& model,
                                                            const std::string& model_path,
                                                            std::optional<int> model_template_index) {
    ModelFileChecker checker(model_path);
    const ModelClassInfo info = checker.resolveModelClass(model);

    if (!info.is_template) {
        return info.class_name;
    }

    if (!model_template_index.has_value()) {
        throw std::runtime_error(
            "MARTY model class '" + info.class_name + "' resolved from mty_model_name='" + model +
            "' is a template class, but no template index was provided."
        );
    }

    return info.class_name + "<" + std::to_string(*model_template_index) + ">";
}

std::string GeneralModelModifier::modelSignature(const std::string& model,
                                                 const std::string& model_path,
                                                 std::optional<int> model_template_index) {
    return "HYPERISO_MARTY_MODEL_SIGNATURE: "
        + resolveModelInstantiation(model, model_path, model_template_index)
        + "; path=" + normalized_path(model_path)
        + "; fnv1a64=" + stable_file_fingerprint(model_path);
}

bool GeneralModelModifier::usesRegPropSplit() const {
    return this->bsm_split_generation
        && (this->wilson == "C9" || this->wilson == "CP9" || this->wilson == "CP10");
}

bool GeneralModelModifier::usesGenericTreeFirst() const {
    return this->tree_first_fallback && !this->usesRegPropSplit();
}

void GeneralModelModifier::replaceWilsonOrderArgument(std::string& line) {
    const auto call = line.find("computeWilsonCoefficients");
    if (call != std::string::npos) {
        const auto open = line.find('(', call);
        if (open != std::string::npos) {
            const std::pair<const char*, const char*> candidates[] = {
                {"mty::Order::TreeLevel", "hyperiso_marty_order"},
                {"mty::Order::OneLoop", "hyperiso_marty_order"},
                {"TreeLevel", "hyperiso_marty_order"},
                {"OneLoop", "hyperiso_marty_order"},
            };
            for (const auto& [needle, replacement] : candidates) {
                const auto pos = line.find(needle, open + 1);
                if (pos != std::string::npos) {
                    line.replace(pos, std::string(needle).size(), replacement);
                    return;
                }
            }
        }
    }

    // Multi-line calls put the order on its own line.
    const auto first = line.find_first_not_of(" \t");
    const auto last = line.find_last_not_of(" \t");
    const std::string trimmed = first == std::string::npos
        ? std::string{}
        : line.substr(first, last - first + 1);
    if (trimmed == "mty::Order::TreeLevel,"
        || trimmed == "mty::Order::OneLoop,"
        || trimmed == "TreeLevel,"
        || trimmed == "OneLoop,") {
        line.replace(first, last - first + 1, "hyperiso_marty_order,");
    }
}

std::string GeneralModelModifier::makeSmFilterHelper() {
    return R"cpp(
namespace {
std::unordered_set<std::string> hyperiso_marty_sm_particle_names() {
    // Do not instantiate an mty::SM_Model inside a generated BSM executable.
    // MARTY keeps global model state internally; constructing a reference SM
    // model while a THDM/SUSY model is active can corrupt that state and has
    // been observed to segfault during fermion embedding.  The split filter
    // only needs a stable name list, so keep it explicit and local.
    return {
        "A", "Z", "W", "Wp", "Wm", "G", "Gp", "Gm", "G0",
        "h",
        "u", "c", "t", "d", "s", "b",
        "u_L", "c_L", "t_L", "d_L", "s_L", "b_L",
        "u_R", "c_R", "t_R", "d_R", "s_R", "b_R",
        "e", "mu", "tau",
        "e_L", "mu_L", "tau_L",
        "e_R", "mu_R", "tau_R",
        "nu_e", "nu_mu", "nu_tau",
        "ve", "vmu", "vtau"
    };
}

std::vector<mty::Particle> hyperiso_marty_non_sm_particles(mty::Model& model) {
    mty::Model::current = &model;
    const auto sm_particle_names = hyperiso_marty_sm_particle_names();

    std::vector<mty::Particle> non_sm_particles;
    for (const auto& particle : model.getParticles()) {
        const std::string name = std::string(particle->getName());
        if (sm_particle_names.find(name) == sm_particle_names.end()) {
            non_sm_particles.push_back(particle);
        }
    }

    return non_sm_particles;
}

void hyperiso_marty_disable_non_sm_particles(mty::FeynOptions& opts, mty::Model& model) {
    auto non_sm_particles = hyperiso_marty_non_sm_particles(model);
    if (!non_sm_particles.empty()) {
        opts.addFilter(mty::filter::disableParticles(non_sm_particles));
    }
}

bool hyperiso_marty_is_non_sm_particle_name(const std::string& name,
                                            const std::unordered_set<std::string>& sm_particle_names) {
    return sm_particle_names.find(name) == sm_particle_names.end();
}

bool hyperiso_marty_has_non_sm_diagram_particle(mty::FeynmanDiagram const& diag,
                                                const std::unordered_set<std::string>& sm_particle_names) {
    for (const auto& particle : diag.getParticles(mty::FeynmanDiagram::DiagramParticleType::Loop)) {
        if (hyperiso_marty_is_non_sm_particle_name(std::string(particle->getName()), sm_particle_names)) {
            return true;
        }
    }
    for (const auto& particle : diag.getParticles(mty::FeynmanDiagram::DiagramParticleType::Mediator)) {
        if (hyperiso_marty_is_non_sm_particle_name(std::string(particle->getName()), sm_particle_names)) {
            return true;
        }
    }
    // MARTY may classify a penguin linker (for example Z_X in b -> s l l)
    // as an External diagram particle even though it is not one of the physical
    // process legs.  Omitting this category silently removes the Z' diagrams.
    for (const auto& particle : diag.getParticles(mty::FeynmanDiagram::DiagramParticleType::External)) {
        if (hyperiso_marty_is_non_sm_particle_name(std::string(particle->getName()), sm_particle_names)) {
            return true;
        }
    }
    return false;
}

void hyperiso_marty_require_non_sm_diagram_particle(mty::FeynOptions& opts) {
    const auto sm_particle_names = hyperiso_marty_sm_particle_names();
    opts.addFilter([sm_particle_names](mty::FeynmanDiagram const& diag) {
        return hyperiso_marty_has_non_sm_diagram_particle(diag, sm_particle_names);
    });
}
} // namespace
)cpp";
}


std::string GeneralModelModifier::makeTreeLevelWilsonHelper() const {
    std::ostringstream configured_order;
    configured_order << "{";
    for (std::size_t i = 0; i < this->tree_fermion_order.size(); ++i) {
        if (i != 0) {
            configured_order << ", ";
        }
        configured_order << this->tree_fermion_order[i];
    }
    configured_order << "}";

    std::string helper = R"cpp(
namespace {
const std::vector<int>& hyperiso_marty_configured_tree_fermion_order()
{
    static const std::vector<int> order = HYPERISO_CONFIGURED_TREE_ORDER;
    return order;
}

std::string hyperiso_marty_fermion_order_label(const std::vector<int>& order)
{
    if (order.empty()) {
        return "template-default";
    }
    std::ostringstream stream;
    for (std::size_t i = 0; i < order.size(); ++i) {
        if (i != 0) {
            stream << '-';
        }
        stream << order[i];
    }
    return stream.str();
}

mty::WilsonSet hyperiso_marty_compute_wilson_coefficients(
    mty::Model& model,
    int order,
    const std::vector<mty::Insertion>& insertions,
    mty::FeynOptions options,
    bool disable_fermion_ordering = false)
{
    int effective_order = order;
#if defined(HYPERISO_MARTY_ORDER_POLICY_TREE_LEVEL_ONLY)
    effective_order = mty::Order::TreeLevel;
#elif defined(HYPERISO_MARTY_ORDER_POLICY_ONE_LOOP_ONLY)
    effective_order = mty::Order::OneLoop;
#endif

    if (effective_order == mty::Order::TreeLevel
        && !hyperiso_marty_configured_tree_fermion_order().empty()) {
        options.setFermionOrder(hyperiso_marty_configured_tree_fermion_order());
        options.orderExternalFermions = false;
    }

    const bool has_explicit_fermion_order = !options.getFermionOrder().empty();
    const bool preserve_explicit_tree_order =
        effective_order == mty::Order::TreeLevel && has_explicit_fermion_order;
    return model.computeWilsonCoefficients(
        effective_order,
        insertions,
        std::move(options),
        disable_fermion_ordering || preserve_explicit_tree_order
    );
}
} // namespace
)cpp";

    const std::string placeholder = "HYPERISO_CONFIGURED_TREE_ORDER";
    helper.replace(helper.find(placeholder), placeholder.size(), configured_order.str());
    return helper;
}

void GeneralModelModifier::replaceWilsonCallWithHelper(std::string& line) {
    const std::string needle = "model.computeWilsonCoefficients(";
    const auto pos = line.find(needle);
    if (pos != std::string::npos) {
        line.replace(pos, needle.size(), "hyperiso_marty_compute_wilson_coefficients(model, ");
    }
}

std::string GeneralModelModifier::orderPolicyPreamble() const {
    switch (this->order_policy) {
    case MartyOrderPolicy::TREE_LEVEL_ONLY:
        return "#define HYPERISO_MARTY_ORDER_POLICY_TREE_LEVEL_ONLY 1\n";
    case MartyOrderPolicy::ONE_LOOP_ONLY:
        return "#define HYPERISO_MARTY_ORDER_POLICY_ONE_LOOP_ONLY 1\n";
    case MartyOrderPolicy::AUTO:
        return "#define HYPERISO_MARTY_ORDER_POLICY_AUTO 1\n";
    }
    return "#define HYPERISO_MARTY_ORDER_POLICY_AUTO 1\n";
}

bool GeneralModelModifier::consumeTreeSafeWilsonCall(std::ofstream& outputFile,
                                                     const std::string& currentLine,
                                                     bool pair_return,
                                                     bool count_graphs) {
    const bool starts_call = currentLine.find("model.computeWilsonCoefficients") != std::string::npos;
    if (!this->buffering_tree_safe_wilson_call && !starts_call) {
        return false;
    }

    if (!this->buffering_tree_safe_wilson_call) {
        this->buffering_tree_safe_wilson_call = true;
        this->tree_safe_wilson_call_lines.clear();
    }
    this->tree_safe_wilson_call_lines.push_back(currentLine);

    if (currentLine.find(");") == std::string::npos) {
        return true;
    }

    emitTreeSafeWilsonCall(outputFile, pair_return, count_graphs);
    this->tree_safe_wilson_call_lines.clear();
    this->buffering_tree_safe_wilson_call = false;
    return true;
}

void GeneralModelModifier::emitTreeSafeWilsonCall(std::ofstream& outputFile,
                                                  bool pair_return,
                                                  bool count_graphs) {
    if (this->tree_safe_wilson_call_lines.empty()) {
        throw std::runtime_error("Cannot emit an empty MARTY Wilson call");
    }

    const std::string& first_line = this->tree_safe_wilson_call_lines.front();
    const auto auto_pos = first_line.find("auto ");
    const auto eq_pos = first_line.find('=', auto_pos == std::string::npos ? 0 : auto_pos + 5);
    const auto call_pos = first_line.find("model.computeWilsonCoefficients", eq_pos);
    if (auto_pos == std::string::npos || eq_pos == std::string::npos
        || call_pos == std::string::npos || eq_pos <= auto_pos + 5) {
        throw std::runtime_error(
            "Cannot rewrite MARTY computeWilsonCoefficients call for a safe tree probe: "
            + first_line
        );
    }

    std::string variable = first_line.substr(auto_pos + 5, eq_pos - (auto_pos + 5));
    const auto variable_first = variable.find_first_not_of(" \t");
    const auto variable_last = variable.find_last_not_of(" \t");
    if (variable_first == std::string::npos || variable_last == std::string::npos) {
        throw std::runtime_error("Cannot resolve MARTY WilsonSet variable from: " + first_line);
    }
    variable = variable.substr(variable_first, variable_last - variable_first + 1);
    const std::string indent = first_line.substr(0, auto_pos);

    auto replace_token = [](std::string& line,
                            const std::string& token,
                            const std::string& replacement) {
        std::size_t pos = 0;
        while ((pos = line.find(token, pos)) != std::string::npos) {
            const bool left_ok = pos == 0
                || !(std::isalnum(static_cast<unsigned char>(line[pos - 1])) || line[pos - 1] == '_');
            const std::size_t right = pos + token.size();
            const bool right_ok = right >= line.size()
                || !(std::isalnum(static_cast<unsigned char>(line[right])) || line[right] == '_');
            if (left_ok && right_ok) {
                line.replace(pos, token.size(), replacement);
                pos += replacement.size();
            } else {
                pos += token.size();
            }
        }
    };

    std::vector<std::string> probe_lines = this->tree_safe_wilson_call_lines;
    {
        std::string& line = probe_lines.front();
        line = indent + "auto hyperiso_marty_tree_probe =" + line.substr(eq_pos + 1);
        const auto method = line.find("computeWilsonCoefficients");
        if (method == std::string::npos) {
            throw std::runtime_error("Cannot locate MARTY Wilson method in buffered call");
        }
        line.replace(method,
                     std::string("computeWilsonCoefficients").size(),
                     "computeAmplitude");
    }
    for (auto& line : probe_lines) {
        replace_token(line, "opts", "hyperiso_marty_tree_options");
    }

    std::vector<std::string> loop_lines = this->tree_safe_wilson_call_lines;
    loop_lines.front() = indent + variable + " =" + loop_lines.front().substr(eq_pos + 1);

    outputFile << indent << "mty::WilsonSet " << variable << ";\n";
    outputFile << indent << "if (hyperiso_marty_order == mty::Order::TreeLevel) {\n";
    outputFile << indent << "    auto hyperiso_marty_tree_options = opts;\n";
    outputFile << indent << "    if (!hyperiso_marty_configured_tree_fermion_order().empty()) {\n";
    outputFile << indent << "        hyperiso_marty_tree_options.setFermionOrder(hyperiso_marty_configured_tree_fermion_order());\n";
    outputFile << indent << "    } else if (!hyperiso_marty_forced_fermion_order.empty()) {\n";
    outputFile << indent << "        hyperiso_marty_tree_options.setFermionOrder(hyperiso_marty_forced_fermion_order);\n";
    outputFile << indent << "    }\n";
    outputFile << indent << "    // Preserve opts.setFermionOrder(...). MARTY's automatic tree ordering can\n";
    outputFile << indent << "    // select a different four-fermion pairing from the one-loop projector.\n";
    outputFile << indent << "    hyperiso_marty_tree_options.orderExternalFermions = false;\n";
    for (const auto& line : probe_lines) {
        outputFile << "    " << line << "\n";
    }
    outputFile << indent << "    if (hyperiso_marty_tree_probe.empty()) {\n";
    if (pair_return) {
        outputFile << indent
                   << "        return std::make_pair(CSL_0, std::size_t{0});\n";
    } else {
        outputFile << indent << "        return CSL_0;\n";
    }
    outputFile << indent << "    }\n";
    outputFile << indent << "    " << variable
               << " = model.getWilsonCoefficients(hyperiso_marty_tree_probe, "
               << "hyperiso_marty_tree_options);\n";
    outputFile << indent << "} else {\n";
    for (const auto& line : loop_lines) {
        outputFile << "    " << line << "\n";
    }
    outputFile << indent << "}\n";
    if (count_graphs) {
        outputFile << indent << "hyperiso_marty_graph_count += "
                   << variable << ".graphs.size();\n";
    }
}

void GeneralModelModifier::modifyLine(std::string& line) {
    if (this->usesRegPropSplit()) {
        // The semileptonic templates contain an order literal as the standalone
        // first argument of computeWilsonCoefficients().  Replace only that
        // argument.  A broad textual replacement also rewrites helper logic such
        // as `order == mty::Order::TreeLevel` into `order == order`, causing all
        // one-loop calls to be treated as tree-level and disabling the C9 linker
        // policy.
        const auto first = line.find_first_not_of(" \t");
        const auto last = line.find_last_not_of(" \t");
        const std::string trimmed = first == std::string::npos
            ? std::string{}
            : line.substr(first, last - first + 1);
        if (trimmed == "mty::Order::TreeLevel,"
            || trimmed == "mty::Order::OneLoop,"
            || trimmed == "TreeLevel,"
            || trimmed == "OneLoop,") {
            line.replace(first, last - first + 1, "hyperiso_marty_order,");
        }
        return;
    }

    if (this->usesGenericTreeFirst()) {
        replaceWilsonOrderArgument(line);
        return;
    }

    replaceWilsonCallWithHelper(line);

    if (line.find("SM_Model sm;") != std::string::npos) {
        line.replace(line.find("SM_Model"), 8, this->model_instantiation);
    }
    else if (line.find("_SM") != std::string::npos) {
        line.replace(line.find("SM"), 2, this->output_model);
    }
}

void GeneralModelModifier::addLine(std::ofstream& outputFile, const std::string& currentLine) {
    if (currentLine.find("dimension6Operator") != std::string::npos) {
        this->tree_four_fermion_projector = true;
    }

    auto is_comment_line = [](const std::string& value) {
        const auto first = value.find_first_not_of(" \t");
        return first != std::string::npos
            && value.compare(first, 2, "//") == 0;
    };

    if (this->skip_old_main) {
        return;
    }

    if (this->usesGenericTreeFirst()) {
        if (currentLine.find("<iostream>") != std::string::npos) {
            outputFile << currentLine << "\n";
            outputFile << "#include <vector>\n";
            outputFile << "#include <utility>\n";
            outputFile << "#include <algorithm>\n";
            outputFile << "#include <sstream>\n";
            if (this->disable_non_sm_particles || this->bsm_split_generation) {
                outputFile << "#include <string>\n";
                outputFile << "#include <unordered_set>\n";
            }
            outputFile << "#include \"" + this->model_path + "\"\n";
            outputFile << "#include \"" + this->marty_path + "\"\n";
            outputFile << "// " << modelSignature(this->target_model, this->model_path, this->model_template_index) << "\n";
            outputFile << orderPolicyPreamble();
            outputFile << "// HYPERISO_MARTY_TREE_FIRST: policy-aware TreeLevel/OneLoop selection\n";
            return;
        }

        if (currentLine.find("using namespace sm_input;") != std::string::npos) {
            outputFile << currentLine << "\n";
            outputFile << makeTreeLevelWilsonHelper() << "\n";
            if (this->disable_non_sm_particles || this->bsm_split_generation) {
                outputFile << makeSmFilterHelper() << "\n";
            }
            return;
        }

        if (currentLine.find("int calculate") != std::string::npos) {
            std::string signature = currentLine;
            const auto int_pos = signature.find("int ");
            const auto open = signature.find('(', int_pos == std::string::npos ? 0 : int_pos + 4);
            const auto close = signature.rfind(')');
            if (int_pos == std::string::npos || open == std::string::npos
                || close == std::string::npos || close <= open) {
                throw std::runtime_error(
                    "Cannot rewrite MARTY calculate signature for tree-first mode: " + currentLine
                );
            }
            this->generic_builder_name = signature.substr(int_pos + 4, open - (int_pos + 4));
            // Insert the new argument before changing the return type. Replacing
            // `int` by the one-character-longer `Expr` first invalidates the
            // previously computed closing-parenthesis index and used to turn
            // `gauge` / `hyperiso_marty_order` into `gaug` /
            // `hyperiso_marty_ordere` in every generic TreeLevel-first wrapper.
            signature.insert(close, ", mty::Order hyperiso_marty_order, const std::vector<int>& hyperiso_marty_forced_fermion_order = {}");
            signature.replace(int_pos, 3, "Expr");
            outputFile << signature << "\n";
            this->inside_calculate_function = true;
            this->expression_returned = false;
            return;
        }

        if (this->inside_calculate_function
            && currentLine.find("FeynOptions opts;") != std::string::npos) {
            outputFile << currentLine << "\n";
            if (this->disable_non_sm_particles) {
                outputFile << "    hyperiso_marty_disable_non_sm_particles(opts, model);\n";
            } else if (this->bsm_split_generation && !this->full_target_generation) {
                outputFile << "    hyperiso_marty_require_non_sm_diagram_particle(opts);\n";
            }
            return;
        }

        if (this->inside_calculate_function
            && consumeTreeSafeWilsonCall(outputFile, currentLine, false, false)) {
            return;
        }

        if (this->inside_calculate_function
            && (currentLine.find("[[maybe_unused]] int sysres") != std::string::npos
                || currentLine.find("mty::Library wilsonLib") != std::string::npos
                || currentLine.find("wilsonLib.cleanExistingSources") != std::string::npos
                || currentLine.find("defineLibPath(wilsonLib)") != std::string::npos
                || currentLine.find("wilsonLib.print") != std::string::npos)) {
            return;
        }

        if (this->inside_calculate_function
            && currentLine.find("wilsonLib.addFunction") != std::string::npos) {
            if (!this->expression_returned) {
                const auto comma = currentLine.find_last_of(',');
                const auto close = currentLine.rfind(')');
                if (comma == std::string::npos || close == std::string::npos || close <= comma) {
                    throw std::runtime_error(
                        "Cannot rewrite MARTY addFunction line for tree-first mode: " + currentLine
                    );
                }
                const std::string expr = currentLine.substr(comma + 1, close - comma - 1);
                outputFile << "    return " << expr << ";\n";
                this->expression_returned = true;
            }
            return;
        }

        if (this->inside_calculate_function
            && currentLine.find("return 0;") != std::string::npos) {
            if (!this->expression_returned) {
                outputFile << "    return CSL_0;\n";
                this->expression_returned = true;
            }
            return;
        }

        if (this->inside_calculate_function && currentLine == "}") {
            this->inside_calculate_function = false;
            this->expression_returned = false;
            outputFile << currentLine << "\n";
            return;
        }

        if (currentLine.find("int main") != std::string::npos) {
            if (this->generic_builder_name.empty()) {
                throw std::runtime_error(
                    "MARTY tree-first generation reached main before calculate function"
                );
            }
            outputFile << "int main() {\n";
            outputFile << "    // MARTY SM input constants carry default numerical values.  They must\n";
            outputFile << "    // be undefined before the target model is constructed; doing this only\n";
            outputFile << "    // inside calculate_*() is too late when a custom model simplifies its\n";
            outputFile << "    // Lagrangian with theta_W, e_em, or other SM inputs in the constructor.\n";
            outputFile << "    undefineNumericalValues();\n";
            outputFile << "    " << this->model_instantiation << " model;\n";
            outputFile << "    std::vector<int> hyperiso_marty_selected_fermion_order = hyperiso_marty_configured_tree_fermion_order();\n";
            outputFile << "    Expr hyperiso_marty_tree = CSL_0;\n";
            if (this->order_policy != MartyOrderPolicy::ONE_LOOP_ONLY) {
                outputFile << "    hyperiso_marty_tree = " << this->generic_builder_name
                           << "(model, gauge::Type::Feynman, mty::Order::TreeLevel, hyperiso_marty_configured_tree_fermion_order());\n";
            }
            if (this->order_policy == MartyOrderPolicy::TREE_LEVEL_ONLY) {
                outputFile << "    Expr hyperiso_marty_selected = hyperiso_marty_tree;\n";
                outputFile << "    const char* hyperiso_marty_selected_order = \"TreeLevelOnly\";\n";
            } else if (this->order_policy == MartyOrderPolicy::ONE_LOOP_ONLY) {
                outputFile << "    Expr hyperiso_marty_selected = " << this->generic_builder_name
                           << "(model, gauge::Type::Feynman, mty::Order::OneLoop, {});\n";
                outputFile << "    const char* hyperiso_marty_selected_order = \"OneLoopOnly\";\n";
            } else {
                outputFile << "    const bool hyperiso_marty_use_tree = hyperiso_marty_tree != CSL_0;\n";
                outputFile << "    Expr hyperiso_marty_selected = hyperiso_marty_tree;\n";
                outputFile << "    if (!hyperiso_marty_use_tree) {\n";
                outputFile << "        hyperiso_marty_selected = " << this->generic_builder_name
                           << "(model, gauge::Type::Feynman, mty::Order::OneLoop, {});\n";
                outputFile << "    }\n";
                outputFile << "    const char* hyperiso_marty_selected_order = hyperiso_marty_use_tree ? \"TreeLevel\" : \"OneLoop\";\n";
            }
            outputFile << "    std::cout << \"[MARTY " << this->wilson
                       << "] selected order=\" << hyperiso_marty_selected_order"
                       << " << \", fermion-order=\" << hyperiso_marty_fermion_order_label(hyperiso_marty_selected_fermion_order)"
                       << " << std::endl;\n";
            outputFile << "    [[maybe_unused]] int sysres = system(\"rm -rf libs/"
                       << this->wilson << "_" << this->output_model << "\");\n";
            outputFile << "    mty::Library wilsonLib(\"" << this->wilson << "_"
                       << this->output_model << "\", \"libs\");\n";
            outputFile << "    wilsonLib.cleanExistingSources();\n";
            outputFile << "    wilsonLib.addFunction(\"" << this->wilson
                       << "\", hyperiso_marty_selected);\n";
            outputFile << "    defineLibPath(wilsonLib);\n";
            outputFile << "    wilsonLib.print();\n";
            outputFile << "    return 0;\n";
            outputFile << "}\n";
            this->skip_old_main = true;
            return;
        }

        outputFile << currentLine << "\n";
        return;
    }

    if (this->usesRegPropSplit()) {
        if (currentLine.find("<iostream>") != std::string::npos) {
            outputFile << currentLine << "\n";
            outputFile << "#include <string>\n";
            outputFile << "#include <unordered_set>\n";
            outputFile << "#include <vector>\n";
            outputFile << "#include <algorithm>\n";
            outputFile << "#include <sstream>\n";
            outputFile << "#include <cstddef>\n";
            outputFile << "#include <utility>\n";
            outputFile << "#include \"" + this->model_path + "\"" << "\n";
            outputFile << "#include \"" + this->marty_path + "\"" << "\n";
            outputFile << "// " << modelSignature(this->target_model, this->model_path, this->model_template_index) << "\n";
            outputFile << orderPolicyPreamble();
            if (this->full_target_generation) {
                outputFile << "// HYPERISO_MARTY_TARGET_SPLIT: complete target-model diagrams in "
                           << this->model_instantiation << "\n";
            } else {
                outputFile << "// HYPERISO_MARTY_BSM_SPLIT: diagrams with at least one non-SM diagram particle in "
                           << this->model_instantiation << "\n";
            }
            outputFile << "// HYPERISO_MARTY_BSM_SPLIT_ABI: model-split-v27\n";
            return;
        }

        if (currentLine.find("using namespace sm_input;") != std::string::npos) {
            outputFile << currentLine << "\n";
            outputFile << makeTreeLevelWilsonHelper() << "\n";
            outputFile << makeSmFilterHelper() << "\n";
            return;
        }

        if (currentLine.find("int calculate") != std::string::npos) {
            outputFile << "std::pair<Expr, std::size_t> hyperiso_marty_build_" << this->wilson
                       << "(Model &model, gauge::Type gauge, mty::Order hyperiso_marty_order, "
                       << "bool hyperiso_marty_sm_like_filter = false, "
                       << "const std::vector<int>& hyperiso_marty_forced_fermion_order = {}) {\n";
            outputFile << "    std::size_t hyperiso_marty_graph_count = 0;\n";
            outputFile << "    hyperiso_marty_set_semileptonic_order(hyperiso_marty_order);\n";
            this->inside_calculate_function = true;
            this->expression_returned = false;
            return;
        }

        if (this->inside_calculate_function && currentLine.find("FeynOptions opts;") != std::string::npos) {
            outputFile << currentLine << "\n";
            outputFile << "    if (hyperiso_marty_sm_like_filter) {\n";
            outputFile << "        hyperiso_marty_disable_non_sm_particles(opts, model);\n";
            if (!this->full_target_generation) {
                outputFile << "    } else {\n";
            }
            if (!this->full_target_generation) {
                outputFile << "        if (hyperiso_marty_order != mty::Order::TreeLevel) {\n";
                outputFile << "            hyperiso_marty_require_non_sm_diagram_particle(opts);\n";
                outputFile << "        }\n";
            }
            outputFile << "    }\n";
            return;
        }

        if (this->inside_calculate_function
            && consumeTreeSafeWilsonCall(outputFile, currentLine, true, true)) {
            return;
        }

        if (is_comment_line(currentLine)) {
            outputFile << currentLine << "\n";
            return;
        }

        if (this->inside_calculate_function && currentLine.find("[[maybe_unused]] int sysres") != std::string::npos) {
            return;
        }
        if (this->inside_calculate_function && currentLine.find("mty::Library wilsonLib") != std::string::npos) {
            return;
        }
        if (this->inside_calculate_function && currentLine.find("wilsonLib.cleanExistingSources") != std::string::npos) {
            return;
        }
        if (this->inside_calculate_function && currentLine.find("defineLibPath(wilsonLib)") != std::string::npos) {
            return;
        }
        if (this->inside_calculate_function && currentLine.find("wilsonLib.print") != std::string::npos) {
            return;
        }
        if (this->inside_calculate_function && currentLine.find("wilsonLib.addFunction") != std::string::npos) {
            if (!this->expression_returned) {
                const auto comma = currentLine.find_last_of(',');
                const auto close = currentLine.rfind(')');
                if (comma == std::string::npos || close == std::string::npos || close <= comma) {
                    throw std::runtime_error("Cannot rewrite MARTY addFunction line for BSM split: " + currentLine);
                }
                std::string expr = currentLine.substr(comma + 1, close - comma - 1);
                outputFile << "    return std::make_pair(" << expr << ", hyperiso_marty_graph_count);\n";
                this->expression_returned = true;
            }
            return;
        }
        if (this->inside_calculate_function && currentLine.find("return 0;") != std::string::npos) {
            if (!this->expression_returned) {
                outputFile << "    return std::make_pair(CSL_0, hyperiso_marty_graph_count);\n";
                this->expression_returned = true;
            }
            return;
        }
        if (this->inside_calculate_function && currentLine == "}") {
            this->inside_calculate_function = false;
            this->expression_returned = false;
            outputFile << currentLine << "\n";
            return;
        }

        if (currentLine.find("int main") != std::string::npos) {
            const bool split_sm_components = false;
            const bool split_linker_components = (this->wilson == "CP10");
            outputFile << "int main() {\n";
            // A fresh C9 model is needed only when AUTO actually falls back
            // from TreeLevel to OneLoop. Constructing both instances before the
            // tree probe changes MARTY's process-global current model and can
            // invalidate an otherwise correct tree projection.
            const bool isolate_c9_auto_fallback =
                this->wilson == "C9" && this->order_policy == MartyOrderPolicy::AUTO;
            const std::string tree_model_name =
                isolate_c9_auto_fallback ? "tree_model" : "model";
            const std::string loop_model_name =
                isolate_c9_auto_fallback ? "loop_model" : "model";
            outputFile << "    // Keep all MARTY SM inputs symbolic from the beginning of target-model\n";
            outputFile << "    // construction.  Otherwise constructor-time simplification can bake in\n";
            outputFile << "    // MARTY defaults before the numerical LHA values are injected.\n";
            outputFile << "    undefineNumericalValues();\n";
            if (isolate_c9_auto_fallback) {
                outputFile << "    " << this->model_instantiation << " tree_model;\n";
            } else {
                outputFile << "    " << this->model_instantiation << " model;\n";
            }
            outputFile << "    std::vector<int> hyperiso_marty_selected_fermion_order = hyperiso_marty_configured_tree_fermion_order();\n";
            outputFile << "    // Semileptonic BSM matching is tree-first.  A non-zero tree-level\n";
            outputFile << "    // coefficient suppresses the one-loop calculation entirely; otherwise\n";
            outputFile << "    // MARTY falls back to the one-loop split-reg_prop path.\n";
            outputFile << "    hyperiso_marty_set_c9_linker_selection(HyperisoMartyC9LinkerSelection::NonPhotonVector);\n";
            outputFile << "#if defined(HYPERISO_MARTY_ORDER_POLICY_ONE_LOOP_ONLY)\n";
            outputFile << "    auto hyperiso_marty_bsm_tree = std::make_pair(CSL_0, std::size_t{0});\n";
            outputFile << "    const bool hyperiso_marty_use_tree_level = false;\n";
            outputFile << "#else\n";
            outputFile << "    auto hyperiso_marty_bsm_tree = hyperiso_marty_build_" << this->wilson
                       << "(" << tree_model_name
                       << ", gauge::Type::Feynman, mty::Order::TreeLevel, false, hyperiso_marty_configured_tree_fermion_order());\n";
            outputFile << "#if defined(HYPERISO_MARTY_ORDER_POLICY_TREE_LEVEL_ONLY)\n";
            outputFile << "    const bool hyperiso_marty_use_tree_level = true;\n";
            outputFile << "#else\n";
            outputFile << "    const bool hyperiso_marty_use_tree_level = hyperiso_marty_bsm_tree.first != CSL_0;\n";
            outputFile << "#endif\n";
            outputFile << "#endif\n";
            outputFile << "    const char* hyperiso_marty_selected_order = hyperiso_marty_use_tree_level ? \"TreeLevel\" : \"OneLoop\";\n";
            outputFile << "    Expr hyperiso_marty_bsm = hyperiso_marty_bsm_tree.first;\n";
            outputFile << "    Expr hyperiso_marty_bsm_photon = CSL_0;\n";
            outputFile << "    Expr hyperiso_marty_bsm_scalar = CSL_0;\n";
            outputFile << "    Expr hyperiso_marty_bsm_vector = hyperiso_marty_bsm_tree.first;\n";
            outputFile << "    std::size_t hyperiso_marty_non_photon_graph_count = hyperiso_marty_bsm_tree.second;\n";
            outputFile << "    std::size_t hyperiso_marty_photon_graph_count = 0;\n";
            outputFile << "    std::size_t hyperiso_marty_scalar_graph_count = 0;\n";
            outputFile << "    std::size_t hyperiso_marty_vector_graph_count = hyperiso_marty_bsm_tree.second;\n";
            outputFile << "    if (!hyperiso_marty_use_tree_level) {\n";
            if (isolate_c9_auto_fallback) {
                outputFile << "        // Instantiate the fallback model only after the tree probe.\n";
                outputFile << "        // Re-assert symbolic SM inputs in case model construction or the\n";
                outputFile << "        // tree probe assigned a process-global MARTY default.\n";
                outputFile << "        undefineNumericalValues();\n";
                outputFile << "        " << this->model_instantiation << " loop_model;\n";
            }
            outputFile << "        hyperiso_marty_set_c9_linker_selection(HyperisoMartyC9LinkerSelection::NonPhotonVector);\n";
            outputFile << "        auto hyperiso_marty_bsm_loop = hyperiso_marty_build_" << this->wilson
                       << "(" << loop_model_name
                       << ", gauge::Type::Feynman, mty::Order::OneLoop, false, {});\n";
            outputFile << "        hyperiso_marty_bsm = hyperiso_marty_bsm_loop.first;\n";
            outputFile << "        hyperiso_marty_non_photon_graph_count = hyperiso_marty_bsm_loop.second;\n";
            outputFile << "        hyperiso_marty_set_c9_linker_selection(HyperisoMartyC9LinkerSelection::PhotonOnly);\n";
            outputFile << "        auto hyperiso_marty_bsm_photon_loop = hyperiso_marty_build_" << this->wilson
                       << "(" << loop_model_name
                       << ", gauge::Type::Feynman, mty::Order::OneLoop, false, {});\n";
            outputFile << "        hyperiso_marty_bsm_photon = hyperiso_marty_bsm_photon_loop.first;\n";
            outputFile << "        hyperiso_marty_photon_graph_count = hyperiso_marty_bsm_photon_loop.second;\n";
            if (split_linker_components) {
                outputFile << "        hyperiso_marty_set_c9_linker_selection(HyperisoMartyC9LinkerSelection::ScalarOnly);\n";
                outputFile << "        auto hyperiso_marty_bsm_scalar_loop = hyperiso_marty_build_" << this->wilson
                           << "(" << loop_model_name
                       << ", gauge::Type::Feynman, mty::Order::OneLoop, false, {});\n";
                outputFile << "        hyperiso_marty_bsm_scalar = hyperiso_marty_bsm_scalar_loop.first;\n";
                outputFile << "        hyperiso_marty_scalar_graph_count = hyperiso_marty_bsm_scalar_loop.second;\n";
                outputFile << "        hyperiso_marty_set_c9_linker_selection(HyperisoMartyC9LinkerSelection::VectorOnly);\n";
                outputFile << "        auto hyperiso_marty_bsm_vector_loop = hyperiso_marty_build_" << this->wilson
                           << "(" << loop_model_name
                       << ", gauge::Type::Feynman, mty::Order::OneLoop, false, {});\n";
                outputFile << "        hyperiso_marty_bsm_vector = hyperiso_marty_bsm_vector_loop.first;\n";
                outputFile << "        hyperiso_marty_vector_graph_count = hyperiso_marty_bsm_vector_loop.second;\n";
            } else {
                outputFile << "        hyperiso_marty_bsm_vector = hyperiso_marty_bsm;\n";
                outputFile << "        hyperiso_marty_vector_graph_count = hyperiso_marty_non_photon_graph_count;\n";
            }
            outputFile << "    }\n";
            if (split_linker_components) {
                outputFile << "    else {\n";
                outputFile << "        // At tree level the penguin linker split is disabled.  Expose the\n";
                outputFile << "        // complete tree coefficient as VECTOR for diagnostics and do not\n";
                outputFile << "        // trigger additional tree or one-loop calculations.\n";
                outputFile << "        hyperiso_marty_bsm_vector = hyperiso_marty_bsm;\n";
                outputFile << "        hyperiso_marty_vector_graph_count = hyperiso_marty_bsm_tree.second;\n";
                outputFile << "    }\n";
            }
            outputFile << "    std::cout << \"[MARTY " << this->wilson << "] "
                       << (this->full_target_generation ? "target" : "BSM")
                       << " selected order=\" << hyperiso_marty_selected_order"
                       << " << \", tree=\" << hyperiso_marty_bsm_tree.second"
                       << " << \", non-photon=\" << hyperiso_marty_non_photon_graph_count"
                       << " << \", photon=\" << hyperiso_marty_photon_graph_count";
            if (split_linker_components) {
                outputFile << " << \", scalar=\" << hyperiso_marty_scalar_graph_count"
                           << " << \", vector=\" << hyperiso_marty_vector_graph_count";
            }
            outputFile << " << \", fermion-order=\" << hyperiso_marty_fermion_order_label(hyperiso_marty_selected_fermion_order)"
                       << " << std::endl;\n";
            if (split_sm_components) {
                outputFile << "    hyperiso_marty_set_c9_linker_selection(HyperisoMartyC9LinkerSelection::NonPhotonVector);\n";
                outputFile << "    auto hyperiso_marty_sm_loop = hyperiso_marty_build_" << this->wilson
                           << "(model, gauge::Type::Feynman, mty::Order::OneLoop, true, {});\n";
                outputFile << "    Expr hyperiso_marty_sm = hyperiso_marty_sm_loop.first;\n";
                outputFile << "    hyperiso_marty_set_c9_linker_selection(HyperisoMartyC9LinkerSelection::PhotonOnly);\n";
                outputFile << "    auto hyperiso_marty_sm_photon_loop = hyperiso_marty_build_" << this->wilson
                           << "(model, gauge::Type::Feynman, mty::Order::OneLoop, true, {});\n";
                outputFile << "    Expr hyperiso_marty_sm_photon = hyperiso_marty_sm_photon_loop.first;\n";
                outputFile << "    hyperiso_marty_set_c9_linker_selection(HyperisoMartyC9LinkerSelection::ScalarOnly);\n";
                outputFile << "    auto hyperiso_marty_sm_scalar_loop = hyperiso_marty_build_" << this->wilson
                           << "(model, gauge::Type::Feynman, mty::Order::OneLoop, true, {});\n";
                outputFile << "    Expr hyperiso_marty_sm_scalar = hyperiso_marty_sm_scalar_loop.first;\n";
                outputFile << "    hyperiso_marty_set_c9_linker_selection(HyperisoMartyC9LinkerSelection::VectorOnly);\n";
                outputFile << "    auto hyperiso_marty_sm_vector_loop = hyperiso_marty_build_" << this->wilson
                           << "(model, gauge::Type::Feynman, mty::Order::OneLoop, true, {});\n";
                outputFile << "    Expr hyperiso_marty_sm_vector = hyperiso_marty_sm_vector_loop.first;\n";
            } else {
                outputFile << "    Expr hyperiso_marty_sm = CSL_0;\n";
                outputFile << "    Expr hyperiso_marty_sm_photon = CSL_0;\n";
                outputFile << "    Expr hyperiso_marty_sm_scalar = CSL_0;\n";
                outputFile << "    Expr hyperiso_marty_sm_vector = CSL_0;\n";
            }
            outputFile << "    hyperiso_marty_set_c9_linker_selection(HyperisoMartyC9LinkerSelection::NonPhotonVector);\n";
            outputFile << "    [[maybe_unused]] int sysres = system(\"rm -rf libs/" << this->wilson << "_" << this->output_model << "\");\n";
            outputFile << "    mty::Library wilsonLib(\"" << this->wilson << "_" << this->output_model << "\", \"libs\");\n";
            outputFile << "    wilsonLib.cleanExistingSources();\n";
            outputFile << "    wilsonLib.addFunction(\"" << this->wilson << "\", hyperiso_marty_bsm);\n";
            outputFile << "    wilsonLib.addFunction(\"" << this->wilson << "_A\", hyperiso_marty_bsm_photon);\n";
            outputFile << "    wilsonLib.addFunction(\"" << this->wilson << "_SCALAR\", hyperiso_marty_bsm_scalar);\n";
            outputFile << "    wilsonLib.addFunction(\"" << this->wilson << "_VECTOR\", hyperiso_marty_bsm_vector);\n";
            outputFile << "    wilsonLib.addFunction(\"" << this->wilson << "_SM\", hyperiso_marty_sm);\n";
            outputFile << "    wilsonLib.addFunction(\"" << this->wilson << "_SM_A\", hyperiso_marty_sm_photon);\n";
            outputFile << "    wilsonLib.addFunction(\"" << this->wilson << "_SM_SCALAR\", hyperiso_marty_sm_scalar);\n";
            outputFile << "    wilsonLib.addFunction(\"" << this->wilson << "_SM_VECTOR\", hyperiso_marty_sm_vector);\n";
            outputFile << "    wilsonLib.addFunction(\"" << this->wilson << "_TOT\", hyperiso_marty_sm + hyperiso_marty_bsm);\n";
            outputFile << "    wilsonLib.addFunction(\"" << this->wilson << "_TOT_A\", hyperiso_marty_sm_photon + hyperiso_marty_bsm_photon);\n";
            outputFile << "    wilsonLib.addFunction(\"" << this->wilson << "_TOT_SCALAR\", hyperiso_marty_sm_scalar + hyperiso_marty_bsm_scalar);\n";
            outputFile << "    wilsonLib.addFunction(\"" << this->wilson << "_TOT_VECTOR\", hyperiso_marty_sm_vector + hyperiso_marty_bsm_vector);\n";
            outputFile << "    defineLibPath(wilsonLib);\n";
            outputFile << "    wilsonLib.print();\n";
            outputFile << "    return 0;\n";
            outputFile << "}\n";
            this->skip_old_main = true;
            return;
        }

        outputFile << currentLine << "\n";
        return;
    }

    if (currentLine.find("<iostream>") != std::string::npos) {
        outputFile << currentLine << "\n";
        outputFile << "#include <vector>\n";
        outputFile << "#include <utility>\n";
        outputFile << "#include <algorithm>\n";
        outputFile << "#include <sstream>\n";
        if (this->disable_non_sm_particles || this->bsm_split_generation) {
            outputFile << "#include <string>\n";
            outputFile << "#include <unordered_set>\n";
        }
        outputFile << "#include \"" + this->model_path + "\"" << "\n";
        outputFile << "#include \"" + this->marty_path + "\"" << "\n";
        outputFile << "// " << modelSignature(this->target_model, this->model_path, this->model_template_index) << "\n";
        outputFile << orderPolicyPreamble();
        if (this->disable_non_sm_particles) {
            outputFile << "// HYPERISO_MARTY_SM_LIKE_FILTER: disable non-SM particles in "
                       << this->model_instantiation << "\n";
        } else if (this->bsm_split_generation) {
            outputFile << "// HYPERISO_MARTY_BSM_ONLY_FILTER: require a non-SM diagram particle in "
                       << this->model_instantiation << "\n";
        }
    }
    else if (currentLine.find("using namespace sm_input;") != std::string::npos) {
        outputFile << currentLine << "\n";
        outputFile << makeTreeLevelWilsonHelper() << "\n";
        if (this->disable_non_sm_particles || this->bsm_split_generation) {
            outputFile << makeSmFilterHelper() << "\n";
        }
    }
    else if (this->disable_non_sm_particles
             && currentLine.find("FeynOptions opts;") != std::string::npos) {
        outputFile << currentLine << "\n";
        outputFile << "    hyperiso_marty_disable_non_sm_particles(opts, model);\n";
    }
    else if (this->bsm_split_generation
             && currentLine.find("FeynOptions opts;") != std::string::npos) {
        outputFile << currentLine << "\n";
        outputFile << "    hyperiso_marty_require_non_sm_diagram_particle(opts);\n";
    }
    else {
        outputFile << currentLine << "\n";
    }

}

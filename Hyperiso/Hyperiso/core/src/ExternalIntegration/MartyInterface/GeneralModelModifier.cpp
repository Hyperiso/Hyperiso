#include "GeneralModelModifier.h"
#include "MartyRuntimeConfig.h"

#include <sstream>
#include <cstddef>
#include <utility>

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
        false
      ) {}

GeneralModelModifier::GeneralModelModifier(std::string wilson,
                                             std::string output_model,
                                             std::string target_model,
                                             std::string model_path,
                                             std::optional<int> model_template_index,
                                             bool disable_non_sm_particles,
                                             bool bsm_split_generation) {
        this->wilson = std::move(wilson);
        this->output_model = std::move(output_model);
        this->target_model = std::move(target_model);
        this->model_path = std::move(model_path);
        this->model_template_index = model_template_index;
        this->disable_non_sm_particles = disable_non_sm_particles;
        this->bsm_split_generation = bsm_split_generation;

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
        + resolveModelInstantiation(model, model_path, model_template_index);
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

bool hyperiso_marty_has_non_sm_internal_particle(mty::FeynmanDiagram const& diag,
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
    return false;
}

void hyperiso_marty_require_non_sm_internal_particle(mty::FeynOptions& opts) {
    const auto sm_particle_names = hyperiso_marty_sm_particle_names();
    opts.addFilter([sm_particle_names](mty::FeynmanDiagram const& diag) {
        return hyperiso_marty_has_non_sm_internal_particle(diag, sm_particle_names);
    });
}
} // namespace
)cpp";
}

void GeneralModelModifier::modifyLine(std::string& line) {
    if (this->bsm_split_generation) {
        auto replace_all = [](std::string& value, const std::string& from, const std::string& to) {
            std::size_t pos = 0;
            while ((pos = value.find(from, pos)) != std::string::npos) {
                value.replace(pos, from.size(), to);
                pos += to.size();
            }
        };
        replace_all(line, "mty::Order::TreeLevel", "hyperiso_marty_order");
        replace_all(line, "mty::Order::OneLoop", "hyperiso_marty_order");
        replace_all(line, "TreeLevel", "hyperiso_marty_order");
        replace_all(line, "OneLoop", "hyperiso_marty_order");
        return;
    }

    if (line.find("SM_Model sm;") != std::string::npos) {
        line.replace(line.find("SM_Model"), 8, this->model_instantiation);
    }
    else if (line.find("_SM") != std::string::npos) {
        line.replace(line.find("SM"), 2, this->output_model);
    }
}

void GeneralModelModifier::addLine(std::ofstream& outputFile, const std::string& currentLine) {
    auto is_comment_line = [](const std::string& value) {
        const auto first = value.find_first_not_of(" \t");
        return first != std::string::npos
            && value.compare(first, 2, "//") == 0;
    };

    if (this->skip_old_main) {
        return;
    }

    if (this->bsm_split_generation) {
        if (currentLine.find("<iostream>") != std::string::npos) {
            outputFile << currentLine << "\n";
            outputFile << "#include <string>\n";
            outputFile << "#include <unordered_set>\n";
            outputFile << "#include <vector>\n";
            outputFile << "#include <cstddef>\n";
            outputFile << "#include <utility>\n";
            outputFile << "#include \"" + this->model_path + "\"" << "\n";
            outputFile << "#include \"" + this->marty_path + "\"" << "\n";
            outputFile << "// " << modelSignature(this->target_model, this->model_path, this->model_template_index) << "\n";
            outputFile << "// HYPERISO_MARTY_BSM_SPLIT: diagrams with at least one non-SM internal particle in "
                       << this->model_instantiation << "\n";
            outputFile << "// HYPERISO_MARTY_BSM_SPLIT_ABI: model-split-v21\n";
            return;
        }

        if (currentLine.find("using namespace sm_input;") != std::string::npos) {
            outputFile << currentLine << "\n";
            outputFile << makeSmFilterHelper() << "\n";
            return;
        }

        if (currentLine.find("int calculate") != std::string::npos) {
            outputFile << "std::pair<Expr, std::size_t> hyperiso_marty_build_" << this->wilson
                       << "(Model &model, gauge::Type gauge, mty::Order hyperiso_marty_order, "
                       << "bool hyperiso_marty_sm_like_filter = false) {\n";
            outputFile << "    std::size_t hyperiso_marty_graph_count = 0;\n";
            this->inside_calculate_function = true;
            this->expression_returned = false;
            this->pending_wilson_graph_count = false;
            this->pending_wilson_set.clear();
            return;
        }

        if (this->inside_calculate_function && currentLine.find("FeynOptions opts;") != std::string::npos) {
            outputFile << currentLine << "\n";
            outputFile << "    if (hyperiso_marty_sm_like_filter) {\n";
            outputFile << "        hyperiso_marty_disable_non_sm_particles(opts, model);\n";
            outputFile << "    } else {\n";
            outputFile << "        hyperiso_marty_require_non_sm_internal_particle(opts);\n";
            outputFile << "    }\n";
            return;
        }

        if (this->inside_calculate_function && this->pending_wilson_graph_count) {
            outputFile << currentLine << "\n";
            if (currentLine.find(");") != std::string::npos) {
                outputFile << "    hyperiso_marty_graph_count += " << this->pending_wilson_set << ".graphs.size();\n";
                this->pending_wilson_graph_count = false;
                this->pending_wilson_set.clear();
            }
            return;
        }

        if (is_comment_line(currentLine)) {
            outputFile << currentLine << "\n";
            return;
        }

        if (this->inside_calculate_function && currentLine.find("model.computeWilsonCoefficients") != std::string::npos) {
            outputFile << currentLine << "\n";
            const auto auto_pos = currentLine.find("auto ");
            const auto eq_pos = currentLine.find('=', auto_pos == std::string::npos ? 0 : auto_pos);
            if (auto_pos != std::string::npos && eq_pos != std::string::npos && eq_pos > auto_pos + 5) {
                this->pending_wilson_set = currentLine.substr(auto_pos + 5, eq_pos - (auto_pos + 5));
                const auto first = this->pending_wilson_set.find_first_not_of(" \t");
                const auto last = this->pending_wilson_set.find_last_not_of(" \t");
                if (first != std::string::npos && last != std::string::npos) {
                    this->pending_wilson_set = this->pending_wilson_set.substr(first, last - first + 1);
                    if (currentLine.find(");") != std::string::npos) {
                        outputFile << "    hyperiso_marty_graph_count += " << this->pending_wilson_set << ".graphs.size();\n";
                    } else {
                        this->pending_wilson_graph_count = true;
                    }
                }
            }
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
            this->pending_wilson_graph_count = false;
            this->pending_wilson_set.clear();
            outputFile << currentLine << "\n";
            return;
        }

        if (currentLine.find("int main") != std::string::npos) {
            const bool split_sm_components = (this->wilson == "CP10");
            outputFile << "int main() {\n";
            outputFile << "    " << this->model_instantiation << " model;\n";
            outputFile << "    hyperiso_marty_set_c9_linker_selection(HyperisoMartyC9LinkerSelection::NonPhotonVector);\n";
            outputFile << "    auto hyperiso_marty_bsm_loop = hyperiso_marty_build_" << this->wilson
                       << "(model, gauge::Type::Feynman, mty::Order::OneLoop, false);\n";
            outputFile << "    Expr hyperiso_marty_bsm = hyperiso_marty_bsm_loop.first;\n";
            outputFile << "    hyperiso_marty_set_c9_linker_selection(HyperisoMartyC9LinkerSelection::PhotonOnly);\n";
            outputFile << "    auto hyperiso_marty_bsm_photon_loop = hyperiso_marty_build_" << this->wilson
                       << "(model, gauge::Type::Feynman, mty::Order::OneLoop, false);\n";
            outputFile << "    Expr hyperiso_marty_bsm_photon = hyperiso_marty_bsm_photon_loop.first;\n";
            if (split_sm_components) {
                outputFile << "    hyperiso_marty_set_c9_linker_selection(HyperisoMartyC9LinkerSelection::NonPhotonVector);\n";
                outputFile << "    auto hyperiso_marty_sm_loop = hyperiso_marty_build_" << this->wilson
                           << "(model, gauge::Type::Feynman, mty::Order::OneLoop, true);\n";
                outputFile << "    Expr hyperiso_marty_sm = hyperiso_marty_sm_loop.first;\n";
                outputFile << "    hyperiso_marty_set_c9_linker_selection(HyperisoMartyC9LinkerSelection::PhotonOnly);\n";
                outputFile << "    auto hyperiso_marty_sm_photon_loop = hyperiso_marty_build_" << this->wilson
                           << "(model, gauge::Type::Feynman, mty::Order::OneLoop, true);\n";
                outputFile << "    Expr hyperiso_marty_sm_photon = hyperiso_marty_sm_photon_loop.first;\n";
            } else {
                outputFile << "    Expr hyperiso_marty_sm = CSL_0;\n";
                outputFile << "    Expr hyperiso_marty_sm_photon = CSL_0;\n";
            }
            outputFile << "    hyperiso_marty_set_c9_linker_selection(HyperisoMartyC9LinkerSelection::NonPhotonVector);\n";
            outputFile << "    [[maybe_unused]] int sysres = system(\"rm -rf libs/" << this->wilson << "_" << this->output_model << "\");\n";
            outputFile << "    mty::Library wilsonLib(\"" << this->wilson << "_" << this->output_model << "\", \"libs\");\n";
            outputFile << "    wilsonLib.cleanExistingSources();\n";
            outputFile << "    wilsonLib.addFunction(\"" << this->wilson << "\", hyperiso_marty_bsm);\n";
            outputFile << "    wilsonLib.addFunction(\"" << this->wilson << "_A\", hyperiso_marty_bsm_photon);\n";
            outputFile << "    wilsonLib.addFunction(\"" << this->wilson << "_SM\", hyperiso_marty_sm);\n";
            outputFile << "    wilsonLib.addFunction(\"" << this->wilson << "_SM_A\", hyperiso_marty_sm_photon);\n";
            outputFile << "    wilsonLib.addFunction(\"" << this->wilson << "_TOT\", hyperiso_marty_sm + hyperiso_marty_bsm);\n";
            outputFile << "    wilsonLib.addFunction(\"" << this->wilson << "_TOT_A\", hyperiso_marty_sm_photon + hyperiso_marty_bsm_photon);\n";
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
        if (this->disable_non_sm_particles) {
            outputFile << "#include <string>\n";
            outputFile << "#include <unordered_set>\n";
            outputFile << "#include <vector>\n";
        }
        outputFile << "#include \"" + this->model_path + "\"" << "\n";
        outputFile << "#include \"" + this->marty_path + "\"" << "\n";
        outputFile << "// " << modelSignature(this->target_model, this->model_path, this->model_template_index) << "\n";
        if (this->disable_non_sm_particles) {
            outputFile << "// HYPERISO_MARTY_SM_LIKE_FILTER: disable non-SM particles in "
                       << this->model_instantiation << "\n";
        }
    }
    else if (this->disable_non_sm_particles
             && currentLine.find("using namespace sm_input;") != std::string::npos) {
        outputFile << currentLine << "\n";
        outputFile << makeSmFilterHelper() << "\n";
    }
    else if (this->disable_non_sm_particles
             && currentLine.find("FeynOptions opts;") != std::string::npos) {
        outputFile << currentLine << "\n";
        outputFile << "    hyperiso_marty_disable_non_sm_particles(opts, model);\n";
    }
    else {
        outputFile << currentLine << "\n";
    }

}

#include "GeneralModelModifier.h"
#include "MartyRuntimeConfig.h"

#include <sstream>

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
        false
      ) {}

GeneralModelModifier::GeneralModelModifier(std::string wilson,
                                             std::string output_model,
                                             std::string target_model,
                                             std::string model_path,
                                             std::optional<int> model_template_index,
                                             bool disable_non_sm_particles) {
        this->wilson = std::move(wilson);
        this->output_model = std::move(output_model);
        this->target_model = std::move(target_model);
        this->model_path = std::move(model_path);
        this->model_template_index = model_template_index;
        this->disable_non_sm_particles = disable_non_sm_particles;

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
std::vector<mty::Particle> hyperiso_marty_non_sm_particles(mty::Model& model) {
    mty::SM_Model sm_reference;
    mty::Model::current = &model;

    std::unordered_set<std::string> sm_particle_names;
    for (const auto& particle : sm_reference.getParticles()) {
        sm_particle_names.emplace(std::string(particle->getName()));
    }

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
} // namespace
)cpp";
}

void GeneralModelModifier::modifyLine(std::string& line) {
    if (line.find("SM_Model sm;") != std::string::npos) {
        line.replace(line.find("SM_Model"), 8, this->model_instantiation);
    }
    else if (line.find("_SM") != std::string::npos) {
        line.replace(line.find("SM"), 2, this->output_model);
    }
}

void GeneralModelModifier::addLine(std::ofstream& outputFile, const std::string& currentLine) {
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

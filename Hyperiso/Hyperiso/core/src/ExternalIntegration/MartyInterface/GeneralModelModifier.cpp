#include "GeneralModelModifier.h"
#include "MartyRuntimeConfig.h"

#include <sstream>

GeneralModelModifier::GeneralModelModifier(std::string wilson,
                                             std::string model,
                                             std::string model_path,
                                             std::optional<int> model_template_index) {
        this->wilson = wilson;
        this->model = std::move(model);
        this->model_path = std::move(model_path);
        this->model_template_index = model_template_index;

        ModelFileChecker checker(this->model_path);
        this->model_class = checker.resolveModelClass(this->model);
        this->model_instantiation = resolveModelInstantiation(
            this->model,
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

void GeneralModelModifier::modifyLine(std::string& line) {
    if (line.find("SM_Model sm;") != std::string::npos) {
        line.replace(line.find("SM_Model"), 8, this->model_instantiation);
    }
    else if (line.find("_SM") != std::string::npos) {
        line.replace(line.find("SM"), 2, this->model);
    }
}

void GeneralModelModifier::addLine(std::ofstream& outputFile, const std::string& currentLine) {
    if (currentLine.find("<iostream>") != std::string::npos) {
        outputFile << currentLine << "\n";
        outputFile << "#include \"" + this->model_path + "\"" << "\n";
        outputFile << "#include \"" + this->marty_path + "\"" << "\n";
        outputFile << "// " << modelSignature(this->model, this->model_path, this->model_template_index) << "\n";
    }
    else {
        outputFile << currentLine << "\n";
    }

}

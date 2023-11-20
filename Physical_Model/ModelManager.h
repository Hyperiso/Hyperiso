#ifndef MODEL_MANAGER_H
#define MODEL_MANAGER_H

#include <memory>
#include <string>
#include "ModelInterface.h"
#include "MSSMModel.h"
#include "NMSSMModel.h"
// Autres modèles

class ModelManager {
public:
    ModelManager();
    void setModel(const std::string& modelName);
    void configureModelParameters(const ModelParameters& params);
    ModelOutput calculateModelPredictions();

private:
    std::unique_ptr<ModelInterface> currentModel;
};

#endif // MODEL_MANAGER_H

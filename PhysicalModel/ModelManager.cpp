// #include "ModelManager.h"

// ModelManager::ModelManager() {}

// void ModelManager::setModel(const std::string& modelName) {
//     if (modelName == "MSSM") {
//         currentModel = std::make_unique<MSSMModel>();
//     } else if (modelName == "NMSSM") {
//         currentModel = std::make_unique<NMSSMModel>();
//     }
//     // Autres modèles
// }

// void ModelManager::configureModelParameters(const ModelParameters& params) {
//     if (currentModel) {
//         currentModel->setParameters(params);
//     }
// }

// ModelOutput ModelManager::calculateModelPredictions() {
//     if (currentModel) {
//         return currentModel->calculatePredictions();
//     }
//     return ModelOutput(); // Retourne un objet vide si aucun modèle n'est défini
// }

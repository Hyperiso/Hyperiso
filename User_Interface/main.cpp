#include <iostream>
#include "ModelManager.h"

int main() {
    ModelManager modelManager;

    std::string modelName;
    std::cout << "Entrez le nom du modèle (MSSM, NMSSM, etc.): ";
    std::cin >> modelName;

    modelManager.setModel(modelName);

    // Choix des paramètres
    ModelParameters params(/* ... */);
    modelManager.configureModelParameters(params);

    // Calculs prédictions
    ModelOutput output = modelManager.calculateModelPredictions();
    std::cout << "Résultats du modèle : " << output.toString() << std::endl;

    return 0;
}

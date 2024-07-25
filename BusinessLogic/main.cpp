#include <iostream>
#include "Chi2Theo.h"

int main() {
    // Charger les paramètres et les observables depuis un fichier JSON
    Chi2Theo& manager = Chi2Theo::getInstance("../DataBase/data_theo.json");

    // Calculer les observables théoriques
    manager.calculate_observables();
    manager.print_observables();
    return 0;

}
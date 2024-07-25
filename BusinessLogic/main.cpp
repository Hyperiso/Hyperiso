#include <iostream>
#include "Chi2Theo.h"
#include "MemoryManager.h"

int main() {
    // Charger les paramètres et les observables depuis un fichier JSON
    auto mm = MemoryManager::GetInstance("Test/testInput.flha", {0, 3});
    mm->init(); 
    Chi2Theo& manager = Chi2Theo::getInstance("../../DataBase/data_theo.json");

    // Calculer les observables théoriques
    // manager.calculate_observables();
    manager.print_observables();
    return 0;

}
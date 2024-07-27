#include <iostream>
#include "Chi2Theo.h"
#include "Chi2Exp.h"
#include "MemoryManager.h"

int main() {
    // Charger les paramètres et les observables depuis un fichier JSON
    auto mm = MemoryManager::GetInstance("Test/testInput.flha", {0, 3});
    mm->init(); 
    Chi2Theo& manager = Chi2Theo::getInstance("../../DataBase/data_theo.json");
    Chi2Exp& manager_exp = Chi2Exp::getInstance("../../DataBase/data_exp.json");
    // Calculer les observables théoriques
    // manager.calculate_observables();
    manager.print_observables();
    manager_exp.print_observables();
    return 0;

}
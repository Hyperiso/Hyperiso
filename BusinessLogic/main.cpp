
#include <iostream>
#include "Chi2Theo.h"
#include "Chi2Exp.h"
#include "MemoryManager.h"
#include "chi2.h"
int main() {
    // Charger les paramètres et les observables depuis un fichier JSON
    auto mm = MemoryManager::GetInstance("Test/testInput.flha", {0, 3});
    mm->init(); 
    Chi2Theo *manager = Chi2Theo::GetInstance("../../DataBase/data_theo.json");
    // Chi2Exp *manager_exp = Chi2Exp::GetInstance("../../DataBase/data_exp.json");

    manager->print_observables();
    // manager_exp->print_observables();
    // manager_exp->print_correlations();
    // manager_exp->print_correlations_matrix();

    Chi2Manager bite(0,0,81., 0);

    // bite.print_inv_cov();

    std::cout << "chi2 : " <<bite.get_chi2() << std::endl;
    return 0;

}
#include "Parametersv2.h"
#include <iostream>
#include <cmath>

int main() {

    auto mm = MemoryManager::GetInstance("Test/InputFiles/testInput.slha", {0,1});  // Initialize program manager with LHA file containing SMINPUTS block
    mm->init();
    auto sm_params = Parameters::GetInstance(0); // SM Model
    auto susy_params = Parameters::GetInstance(1); // SUSY Model
    
    try {
        
        double alpha_s_MZ = pow((*sm_params)("GAUGE", 1), 2) / (4 * M_PI);
        std::cout << "Alpha_s(MZ): " << alpha_s_MZ << std::endl;

        std::cout << "mass up quark : " << (*sm_params)("MASS", 2) << std::endl;
        std::cout << "first element of ckm (real part) : " << (*sm_params)("RECKM", 0) << std::endl;
         double susy_mass = (*susy_params)("MASS", 1000021); // Example PDG code for a SUSY particle
        std::cout << "SUSY Particle Mass: " << susy_mass << std::endl;


    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    try {
        std::cout << (*susy_params)(std::string("GAUGE"), 1) << std::endl;
    } catch(const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
 
    // auto mm2 = MemoryManager::GetInstance("Test/InputFiles/testinput_thdm.lha", {0,2});
    // mm2->init();
    // auto thdm_params = Parameters::GetInstance(2); // THDM Model
    // std::cout << "THDM matrice yu : " << (*thdm_params)("YU", 22) << std::endl;

    return 0;
}
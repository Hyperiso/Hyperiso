#include <iostream>
#include <fstream>
#include "./Physical_Model/Wilson.h"

int main() {
    std::ofstream file("C7_running.csv");
    file << "Q,C7_real,C7_imag\n"; 

    double Q_initial = 1.0;
    double Q_final = 500.0;
    double Q_step = 1.0;

    
    auto strat = std::make_shared<SM_LO_Strategy>();
    WilsonManager* wm = WilsonManager::GetInstance(100.0, strat);

    
    for (double Q = Q_initial; Q <= Q_final; Q += Q_step) {
        wm->setScale(Q);
        complex_t C7 = wm->get(WilsonCoefficient::C7, 0); 

        file << Q << "," << C7.real() << "," << C7.imag() << "\n";
    }

    
    file.close();
    // delete wm; 

    std::cout << "Done writing C7_running.csv" << std::endl;

    return 0;
}
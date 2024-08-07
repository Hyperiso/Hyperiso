#include "QCDParameters.h"
#include <fstream>
#include <iostream>

int main() {

    double Qinit = 1.0; 
    double Qfin = 500.0;
    double step = 1.0;
    double mtop = 172.9;
    double mbot = 4.18; 
    double mcharm = 1.275; 

    QCDParameters qcdParams = QCDParameters();


    std::ofstream file("quark_mass_running.csv");
    file << "Q,RunningMassTop,RunningMassBottom,RunningMassCharm\n";

    for (double Q = Qinit; Q <= Qfin; Q += step) {
        double runningMassTop = qcdParams.running_mass(mtop, Qinit, Q); 
        double runningMassBottom = qcdParams.running_mass(mbot, Qinit, Q);
        double runningMassCharm = qcdParams.running_mass(mcharm, Qinit, Q); 
        file << Q << "," << runningMassTop << "," << runningMassBottom << "," << runningMassCharm << "\n";
    }

    file.close();

    std::cout << "Running de la masse enregistré dans 'quark_mass_running.csv'" << std::endl;

    return 0;
}
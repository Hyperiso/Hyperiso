#include "QCDParameters.h"
#include <fstream>
#include <iostream>

int main() {
    // Initialisation des paramètres
    double Qinit = 1.0; // Échelle d'énergie initiale en GeV
    double Qfin = 500.0; // Échelle d'énergie finale en GeV
    double step = 1.0; // Pas d'échelle d'énergie pour le calcul
    double mtop = 172.9; // Masse du quark top en GeV
    double mbot = 4.18; // Masse du quark bottom en GeV
    double mcharm = 1.275; // Masse du quark charm en GeV

    QCDParameters qcdParams = QCDParameters(0.1, 91, 4.19000000e+00, 1.72900000e+02, 4.19000000e+00, 1.72900000e+02);

    // Ouverture du fichier CSV pour écrire
    std::ofstream file("quark_mass_running.csv");
    file << "Q,RunningMassTop,RunningMassBottom,RunningMassCharm\n"; // En-tête du fichier

    // Boucle sur les échelles d'énergie
    for (double Q = Qinit; Q <= Qfin; Q += step) {
        double runningMassTop = qcdParams.running_mass(mtop, Qinit, Q, mtop, mbot); // Assumer que la fonction accepte ces paramètres
        double runningMassBottom = qcdParams.running_mass(mbot, Qinit, Q, mtop, mbot); // Assumer que la fonction accepte ces paramètres
        double runningMassCharm = qcdParams.running_mass(mcharm, Qinit, Q, mtop, mbot); // Assumer que la fonction accepte ces paramètres
        file << Q << "," << runningMassTop << "," << runningMassBottom << "," << runningMassCharm << "\n";
    }

    // Fermeture du fichier
    file.close();

    std::cout << "Running de la masse enregistré dans 'quark_mass_running.csv'" << std::endl;

    return 0;
}
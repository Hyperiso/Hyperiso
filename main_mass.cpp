#include "QCDParameters.h"
#include <fstream>
#include <iostream>


int main() {
    // Initialisation des paramètres
    double Qinit = 1.0; // Échelle d'énergie initiale en GeV
    double Qfin = 500.0; // Échelle d'énergie finale en GeV
    double step = 1.0; // Pas d'échelle d'énergie pour le calcul
    double mtop = 172.9; // Masse du quark top
    double mbot = 4.18; // Masse du quark bottom
    double quark_mass = mbot; // Exemple : calculer le running pour la masse du quark bottom

    QCDParameters qcdParams = QCDParameters(0.1, 91, 4.19000000e+00, 1.72900000e+02,4.19000000e+00,1.72900000e+02);

    // Ouverture du fichier CSV pour écrire
    std::ofstream file("quark_mass_running.csv");
    file << "Q,RunningMass\n"; // En-tête du fichier

    // Boucle sur les échelles d'énergie
    for (double Q = Qinit; Q <= Qfin; Q += step) {
        double runningMass = qcdParams.running_mass(quark_mass, Qinit, Q, mtop, mbot);
        file << Q << "," << runningMass << "\n";
    }

    // Fermeture du fichier
    file.close();

    std::cout << "Running de la masse enregistré dans 'quark_mass_running.csv'" << std::endl;

    return 0;
}
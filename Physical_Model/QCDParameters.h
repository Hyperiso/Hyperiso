#include <cmath>

class QCDParameters {
public:
    // Paramètres à initialiser
    double mass_Z;        // MZ
    double alphas_MZ;     // alphas(MZ)
    double Lambda3, Lambda4, Lambda5, Lambda6;
    double alphasMZ_Lambda3, alphasMZ_Lambda4, alphasMZ_Lambda5, alphasMZ_Lambda6;
    double mass_c, mass_b, mass_t;
    const double pi = 3.14159265358979323846;

    QCDParameters();

    double alphasRunning(double Q, double Lambda, int nf);
    double DichotomieLambda(double alpha_running, double Q, int nf);
    double runningAlphasCalculation(double Q, double Lambda, int nf);
    double alphas_running(double Q, double mtop, double mbot);
    double computeAlphas(double Q, double Lambda, int nf);

};

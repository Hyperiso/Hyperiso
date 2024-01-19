#include <cmath>

class QCDParameters {
public:
    // Paramètres à initialiser
    int nf {5};
    double mass_Z;        // MZ
    double alphas_MZ;     // alphas(MZ)
    double Lambda5 {0.};
    double Lambda3, Lambda4, Lambda6;
    double alphasMZ_Lambda3, alphasMZ_Lambda4, alphasMZ_Lambda5, alphasMZ_Lambda6;
    double mass_c;
    double mass_b {4.19999981};
    double mass_t {172.399994};
    const double pi {3.14159265358979323846};

    
    QCDParameters();

    double alphasRunning(double Q, double Lambda, int nf);
    double DichotomieLambda(double alpha_running, double Q, int nf);
    double runningAlphasCalculation(double Q);

};

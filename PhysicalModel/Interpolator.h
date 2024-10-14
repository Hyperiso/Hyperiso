#include <complex>

class Interpolator {
public:
    static std::complex<double> linearInterpolation(double Q1, double Q2, double Q_target, 
                                                    std::complex<double> C1, std::complex<double> C2) {
        double t = (Q_target - Q1) / (Q2 - Q1);
        return C1 + t * (C2 - C1);
    }
};
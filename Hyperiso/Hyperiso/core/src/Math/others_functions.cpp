#include <cmath>
#include "Math.h"

double expcor(double x) {
    if (x > 100.0) {
        return std::exp(100.0);
    } else if (x < -100.0) {
        return 0.0;
    } else {
        return std::exp(x);
    }
}

// Kronecker delta
double kron(int x, int y) {
    return (x == y) ? 1.0 : 0.0;
}
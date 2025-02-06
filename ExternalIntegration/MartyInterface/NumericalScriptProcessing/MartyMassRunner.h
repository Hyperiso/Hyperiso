#include "Parameters.h"

class MartyMassRunner {
    double Q_match;
public:
    MartyMassRunner(double Q_match) : Q_match(Q_match) {}

    double get_mt() {
        return QCDHelper::msbar_mass(6, Q_match);
    }

    double get_mb() {
        return QCDHelper::msbar_mass(5, Q_match);
    }

    double get_mc() {
        return QCDHelper::msbar_mass(4, Q_match);
    }
};
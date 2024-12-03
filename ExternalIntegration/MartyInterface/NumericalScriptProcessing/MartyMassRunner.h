#include "Parameters.h"

class MartyMassRunner {
    double Q_match;
public:
    MartyMassRunner(double Q_match) : Q_match(Q_match) {}

    double get_mt() {
        double m_t = (*Parameters::GetInstance(ParameterType::SM))("MASS", 6);
        return Parameters::GetInstance(ParameterType::SM)->running_mass(m_t, m_t, Q_match);
    }

    double get_mb() {
        double m_b = (*Parameters::GetInstance(ParameterType::SM))("MASS", 5);
        return Parameters::GetInstance(ParameterType::SM)->running_mass(m_b, m_b, Q_match);
    }

    double get_mc() {
        double m_c = (*Parameters::GetInstance(ParameterType::SM))("MASS", 4);
        return Parameters::GetInstance(ParameterType::SM)->running_mass(m_c, m_c, Q_match);
    }
};
// EpsilonCalculator.h
#ifndef EPSILONCALCULATOR_H
#define EPSILONCALCULATOR_H
#include "Parameters.h"
#include <cmath>
#include <vector>
#include <map>

class EpsilonCalculator {
protected:
    Parameters* sm = Parameters::GetInstance();
    Parameters* susy = Parameters::GetInstance(1);

    double mu_Q = (*susy)("HMIX",1);
    std::map<int,int> neutralino = {{1, 1000022},{2, 1000023},{3, 1000025},{4, 1000035}};
    static EpsilonCalculator* instance;
    EpsilonCalculator();

public:
    // explicit EpsilonCalculator();
    
    double epsilon_0();
    double epsilon_2() const;
    double epsilon_b();
    double epsilon_bp();
    double epsilon_0p();
    double epsilon_1p() const;

    EpsilonCalculator(const EpsilonCalculator&) = delete;
    void operator=(const EpsilonCalculator&) = delete;

    // Méthode statique pour accéder à l'instance
    static EpsilonCalculator* GetInstance() {
        if (!EpsilonCalculator::instance) {
            EpsilonCalculator::instance = new EpsilonCalculator();
        }
        return EpsilonCalculator::instance;
    }
};


#endif

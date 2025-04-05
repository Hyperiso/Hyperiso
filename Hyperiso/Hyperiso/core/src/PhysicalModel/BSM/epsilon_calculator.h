// EpsilonCalculator.h
#ifndef EPSILONCALCULATOR_H
#define EPSILONCALCULATOR_H
#include "Parameters.h"
#include "QCDHelper.h"
#include "WilsonParamComposer.h"
#include <cmath>
#include <vector>
#include <map>

class EpsilonCalculator {
protected:
    std::shared_ptr<Parameters> sm = Parameters::GetInstance(ParameterType::SM);
    std::shared_ptr<Parameters> susy = Parameters::GetInstance(ParameterType::BSM);

    double mu_Q = (*susy)("HMIX",1);
    
    
    static inline bool initialized {false};
    
    public:
    
    // double epsilon_0();
    // double epsilon_2() const;
    // double epsilon_b();
    // double epsilon_bp();
    // double epsilon_0p();
    // double epsilon_1p() const;

    static void init();

    static void init_epsilon_block();


    static inline WilsonParamComposer composer = WilsonParamComposer();
};


#endif

#ifndef RARE_DECAY_CALCULATOR_H
#define RARE_DECAY_CALCULATOR_H

#include "ModelParameters.h"

class RareDecayCalculator {
public:
    RareDecayCalculator(ModelParameters params);
    double calculateBranchingRatio();
    void updateModelParameters(ModelParameters newParams);

private:
    ModelParameters modelParams;
};

#endif // RARE_DECAY_CALCULATOR_H

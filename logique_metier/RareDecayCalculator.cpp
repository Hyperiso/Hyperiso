#include "RareDecayCalculator.h"

RareDecayCalculator::RareDecayCalculator(ModelParameters params)
    : modelParams(params) {}

double RareDecayCalculator::calculateBranchingRatio() {

    double mass = modelParams.getParticleMass();
    double couplingConstant = modelParams.getCouplingConstant();


    double branchingRatio = mass * couplingConstant / 1000.0;
    return branchingRatio;
}

void RareDecayCalculator::updateModelParameters(ModelParameters newParams) {
    modelParams = newParams;
}

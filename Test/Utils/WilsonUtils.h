#ifndef WILSON_UTILS_H
#define WILSON_UTILS_H

#include <string>
#include <memory>
#include "Wilson.h"
#include "Wilson_susy.h"
#include "Wilson_THDM.h"

void writeCoefficientsToFile(const std::string& strat_name, const std::string& fileName, const std::shared_ptr<InitializationStrategy>& strategy, double Q_match, const std::string& model);
void writeCoefficientsPrimeCQToFile(const std::string& strat_name, const std::string& fileName, const std::shared_ptr<InitializationStrategy>& strategy, double Q_match, const std::string& model);
void runTest(const std::string& strategyName, const std::shared_ptr<InitializationStrategy>& strategy, const std::string& testFile, const std::string& referenceFile,const std::string& model, double tolerance);
#endif // WILSON_UTILS_H

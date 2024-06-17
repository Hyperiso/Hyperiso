#ifndef WILSON_UTILS_H
#define WILSON_UTILS_H

#include <string>
#include <memory>
#include "Wilson.h"

void writeCoefficientsToFile(const std::string& strat_name, const std::string& fileName, const std::shared_ptr<InitializationStrategy>& strategy, double Q_match, const std::string& model);
void writeCoefficientsPrimeCQToFile(const std::string& strat_name, const std::string& fileName, const std::shared_ptr<InitializationStrategy>& strategy, double Q_match, const std::string& model);

#endif // WILSON_UTILS_H

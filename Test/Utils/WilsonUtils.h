#ifndef WILSON_UTILS_H
#define WILSON_UTILS_H

#include <string>
#include <memory>
#include "WilsonManager.h"
#include "Wilson_THDMv2.h"
#include "Wilson_susyv2.h"
// #include "Wilson.h"
// #include "Wilson_susy.h"
// #include "Wilson_THDM.h"

void writeCoefficientsToFile(const std::string& strat_name, const std::string& fileName, double Q_match, const std::string& model);
void writeCoefficientsPrimeCQToFile(const std::string& strat_name, const std::string& fileName, double Q_match, const std::string& model);
void runTest(const std::string& strategyName, const std::string& testFile, const std::string& referenceFile,const std::string& model, double tolerance, bool PrimeCQ = false);

#endif // WILSON_UTILS_H


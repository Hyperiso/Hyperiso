#ifndef TWO_HDMC_CALCULATOR_H
#define TWO_HDMC_CALCULATOR_H

#include "Interface.h"

/**
 * Concrete spectrum calculator backed by the vendored 2HDMC sources.
 *
 * 2HDMC is part of Hyperiso's ExternalIntegration module and is therefore
 * always available in the standard build.
 */
class TwoHDMCalculator : public ICalculator {
public:
    void calculateSpectrum(const std::string& inputFilePath, const std::string& outputFilePath) override;
};

#endif // TWO_HDMC_CALCULATOR_H

#ifdef BUILD_WITH_SOFTSUSY
#ifndef SOFT_SUSY_H
#define SOFT_SUSY_H

#include <string>
#include <functional> 
#include <unordered_map> 

#include "Interface.h"
#include "config.hpp"

/**
 * @file SoftSusy.h
 * @brief Declares the SoftsusyCalculator concrete spectrum calculator.
 *
 * This header defines ::SoftsusyCalculator, a concrete implementation of
 * ::ICalculator that wraps the SOFTSUSY executable (or library) to
 * compute SUSY spectra from an SLHA / Les Houches input.
 *
 * @ingroup SpectrumCalculationModule
 */

/**
 * @class SoftsusyCalculator
 * @ingroup SpectrumCalculationModule
 * @brief Concrete implementation of ICalculator for SOFTSUSY.
 *
 * SoftsusyCalculator builds and executes a SOFTSUSY command line, taking
 * into account the project’s third-party and assets root paths. It is
 * typically instantiated through ::GeneralCalculatorFactory when
 * ::CalculatorType::Softsusy is requested.
 */
class SoftsusyCalculator : public ICalculator {
public:
    /**
     * @brief Calculates the particle spectrum using SOFTSUSY.
     *
     * The method:
     *  - builds a shell command pointing to `softpoint.x`,
     *  - feeds @p inputFilePath as stdin (possibly prefixed by assets root),
     *  - redirects stdout to @p outputFilePath,
     *  - checks the exit code and logs an error on failure.
     *
     * @param inputFilePath   Path to the input file containing model parameters.
     * @param outputFilePath  Path where the calculated spectrum should be written.
     */
    void calculateSpectrum(const std::string& inputFilePath, const std::string& outputFilePath) override;
};

#endif
#endif

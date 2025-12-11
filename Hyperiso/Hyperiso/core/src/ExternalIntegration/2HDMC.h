#ifdef BUILD_WITH_2HDMC
#ifndef T2HDMC_H
#define T2HDMC_H

#include <unordered_map>
#include <functional>
#include <iostream>

#include "Interface.h"

/**
 * @file 2HDMC.h
 * @brief Declares the TwoHDMCalculator concrete spectrum calculator.
 *
 * This header provides the ::TwoHDMCalculator class, a concrete implementation
 * of ::ICalculator that uses the 2HDMC library to compute spectra in
 * Two-Higgs-Doublet Models (THDM).
 *
 * @ingroup SpectrumCalculationModule
 */

 /**
 * @class TwoHDMCalculator
 * @ingroup SpectrumCalculationModule
 * @brief Concrete implementation of ICalculator using 2HDMC.
 *
 * TwoHDMCalculator:/**
 * @file SoftSusy.h
 * @brief Declares the SoftsusyCalculator concrete spectrum calculator.
 *
 * This header defines ::SoftsusyCalculator, a concrete implementation of
 * ::ICalculator that wraps the SOFTSUSY executable (or library) to
 * compute SUSY spectra from an SLHA / Les Houches input.
 *
 * @ingroup SpectrumCalculationModule
 */
 *  - reads an LHA-like input file describing a THDM point,
 *  - delegates to 2HDMC (via ::THDM) to interpret the file,
 *  - writes out a spectrum / SLHA-style output.
 *
 * It is created by ::GeneralCalculatorFactory when ::CalculatorType::TwoHDM
 * is requested.
 */
class TwoHDMCalculator : public ICalculator {
public:
    /**
     * @brief Calculates the spectrum using 2HDMC.
     *
     * The method:
     *  - builds a ::THDM model instance,
     *  - calls THDM::read_LesHouches() on @p inputFilePath,
     *  - writes the resulting SLHA/LHA output to @p outputFilePath.
     *
     * @param inputFilePath   Path to the input file containing THDM parameters.
     * @param outputFilePath  Path where the calculated spectrum will be written.
     */
    void calculateSpectrum(const std::string& inputFilePath, const std::string& outputFilePath) override;
};

#endif
#endif
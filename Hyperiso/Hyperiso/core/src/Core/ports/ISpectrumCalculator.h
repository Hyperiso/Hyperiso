#ifndef ISPECTRUMCALCULATOR_H
#define ISPECTRUMCALCULATOR_H

#include <filesystem>

#include "Include.h"

namespace fs = std::filesystem;

/**
 * @file ISpectrumCalculator.h
 * @brief Declares the interface for spectrum calculators.
 *
 * This header defines the abstract interface ::ISpectrumCalculator, which
 * encapsulates the notion of a “spectrum calculator” driven by an input
 * LHA file and a chosen physics model.
 *
 * Typical responsibilities:
 *  - read an input LHA file describing model parameters,
 *  - run an external spectrum generator,
 *  - write the resulting spectrum to an output LHA file.
 *
 * @example spectrum_calculator_example.cpp
 * @brief Example usage of SpectrumCalculator to calculate a physics spectrum.
 *
 * @defgroup SpectrumCalculationModule Spectrum Calculation System
 * @brief Provides interfaces and classes to calculate particle physics spectra.
 *
 * ## Related Classes
 * - @ref ISpectrumCalculator
 * - @ref SpectrumCalculator
 */

/**
 * @class ISpectrumCalculator
 * @ingroup SpectrumCalculationModule
 * @brief Interface for spectrum calculators based on input LHA files.
 *
 * Implementations of this interface wrap concrete spectrum generators
 * (e.g. TwoHDM-specific tools, SUSY codes like Softsusy, etc.).
 *
 * The primary entry point is ::calculate_spectrum(), which:
 *  - takes an existing LHA file as input,
 *  - performs spectrum calculation according to the provided Model,
 *  - writes the resulting spectrum to a separate LHA file.
 */
class ISpectrumCalculator {
public:
    /// Virtual destructor for interface.
    virtual ~ISpectrumCalculator() = default;

    /**
     * @brief Calculates the spectrum based on an input LHA file and outputs the result.
     *
     * Implementations are expected to:
     *  - read the input file @p in_lha_path,
     *  - run an external spectrum generator consistent with @p model,
     *  - write the resulting spectrum to @p out_spectrum_path.
     *
     * Existing files at @p out_spectrum_path may be overwritten depending
     * on the concrete implementation.
     *
     * @param in_lha_path         Path to the input LHA file.
     * @param out_spectrum_path   Path where the calculated spectrum should be saved.
     * @param model               The physics model to use for the spectrum calculation.
     */
    virtual void calculate_spectrum(fs::path in_lha_path, fs::path out_spectrum_path, Model model) = 0;

    
};

#endif // ISPECTRUMCALCULATOR_H
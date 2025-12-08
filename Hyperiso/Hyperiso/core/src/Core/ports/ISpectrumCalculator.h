#ifndef ISPECTRUMCALCULATOR_H
#define ISPECTRUMCALCULATOR_H

#include <filesystem>

#include "Include.h"

/**
 * @example spectrum_calculator_example.cpp
 * @brief Example usage of SpectrumCalculator to calculate a physics spectrum.
 * @defgroup SpectrumCalculationModule Spectrum Calculation System
 * @brief Provides interfaces and classes to calculate particle physics spectra.
 *
 * ## Related Classes
 * - @ref ISpectrumCalculator
 * - @ref SpectrumCalculator
 */

namespace fs = std::filesystem;

/**
 * @class ISpectrumCalculator
 * @ingroup SpectrumCalculationModule
 * @brief Interface for spectrum calculators based on input LHA files.
 *
 * Defines the contract for classes that perform spectrum calculation
 * based on different physics models.
 */
class ISpectrumCalculator {
public:

    /**
     * @brief Calculates the spectrum based on an input LHA file and outputs the result.
     * @param in_lha_path Path to the input LHA file.
     * @param out_spectrum_path Path where the calculated spectrum should be saved.
     * @param model The physics model to use for the spectrum calculation.
     */
    virtual void calculate_spectrum(fs::path in_lha_path, fs::path out_spectrum_path, Model model) = 0;

    virtual ~ISpectrumCalculator() = default;
};

#endif // ISPECTRUMCALCULATOR_H
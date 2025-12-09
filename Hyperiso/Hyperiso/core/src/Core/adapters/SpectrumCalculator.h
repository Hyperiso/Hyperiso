#ifndef SPECTRUMCALCULATOR_H
#define SPECTRUMCALCULATOR_H

#include "ISpectrumCalculator.h"
#include "Interface.h"

/**
 * @file SpectrumCalculator.h
 * @brief Concrete spectrum calculator using external tools.
 *
 * This header declares the ::SpectrumCalculator class, which provides a
 * concrete implementation of ::ISpectrumCalculator based on the generic
 * calculator interface defined in `Interface.h`.
 *
 * Depending on the model, it chooses an appropriate CalculatorType and
 * delegates the actual computation to GeneralCalculatorFactory.
 */

/**
 * @class SpectrumCalculator
 * @ingroup SpectrumCalculationModule
 * @brief Concrete implementation of ISpectrumCalculator.
 *
 * SpectrumCalculator wraps calls to external spectrum generators through
 * ::GeneralCalculatorFactory. The mapping between Model and CalculatorType
 * is handled internally:
 *  - Model::THDM  → CalculatorType::TwoHDM
 *  - other BSM    → CalculatorType::Softsusy (current default)
 *
 * It is designed to be plugged into MemoryManager via ISpectrumCalculator,
 * so that spectrum computation can be transparently toggled or replaced.
 */
class SpectrumCalculator : public ISpectrumCalculator {
public:
    /**
     * @brief Calculates the physics spectrum by calling the appropriate external calculator.
     *
     * For a given @p model, the method:
     *  - selects a CalculatorType,
     *  - invokes GeneralCalculatorFactory::executeCommand("calculateSpectrum", ...),
     *  - logs the beginning and successful end of the spectrum computation.
     *
     * @param in_lha_path         Path to the input LHA file.
     * @param out_spectrum_path   Path where the output spectrum should be saved.
     * @param model               The model to be used for the calculation (e.g., THDM, SUSY).
     */
    void calculate_spectrum(fs::path in_lha_path, fs::path out_spectrum_path, Model model) override;
};


#endif // SPECTRUMCALCULATOR_H

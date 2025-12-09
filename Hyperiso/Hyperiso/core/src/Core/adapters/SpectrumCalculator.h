#ifndef SPECTRUMCALCULATOR_H
#define SPECTRUMCALCULATOR_H

#include "ISpectrumCalculator.h"
#include "Interface.h"

/**
 * @class SpectrumCalculator
 * @ingroup SpectrumCalculationModule
 * @brief Concrete implementation of ISpectrumCalculator.
 *
 * Uses an external tool to calculate the physics spectrum based on the specified model.
 */
class SpectrumCalculator : public ISpectrumCalculator {
public:

    /**
     * @brief Calculates the physics spectrum by calling the appropriate external calculator.
     *
     * @param in_lha_path         Path to the input LHA file.
     * @param out_spectrum_path   Path where the output spectrum should be saved.
     * @param model               The model to be used for the calculation (e.g., THDM, SUSY).
     */
    void calculate_spectrum(fs::path in_lha_path, fs::path out_spectrum_path, Model model) override;
};


#endif // SPECTRUMCALCULATOR_H

#include "SpectrumCalculator.h"

void SpectrumCalculator::calculate_spectrum(fs::path in_lha_path, fs::path out_spectrum_path, Model model) {
    LOG_DEBUG("Starting spectrum calculation...");
    CalculatorType calculatorType = model == Model::THDM ? CalculatorType::TwoHDM : CalculatorType::Softsusy;
    GeneralCalculatorFactory::executeCommand(calculatorType, "calculateSpectrum", in_lha_path.string(), out_spectrum_path.string());
    LOG_DEBUG("Spectrum calculation completed successfully");
}
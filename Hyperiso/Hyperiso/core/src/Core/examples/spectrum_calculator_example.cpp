/**
 * @example spectrum_calculator_example.cpp
 * @brief Example usage of SpectrumCalculator to calculate a physics spectrum.
 *
 * @code
 * #include "SpectrumCalculator.h"
 * #include <iostream>
 *
 * int main() {
 *     SpectrumCalculator calc;
 *
 *     // Input and output file paths
 *     fs::path input_lha = "input.lha";     // Path to your input LHA file
 *     fs::path output_spectrum = "spectrum_output.lha"; // Where to save the result
 *
 *     // Specify the model to use for the calculation
 *     Model model = Model::THDM; // or Model::SUSY
 *
 *     // Perform the calculation
 *     calc.calculate_spectrum(input_lha, output_spectrum, model);
 *
 *     std::cout << "Spectrum calculation completed successfully." << std::endl;
 *
 *     return 0;
 * }
 * @endcode
 */

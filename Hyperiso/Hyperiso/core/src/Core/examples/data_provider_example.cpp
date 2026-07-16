/**
 * @example data_provider_example.cpp
 * @brief Example usage of ParameterProvider, CorrelationProvider, and QCDProvider.
 *
 * @code
 * #include "ParameterProvider.h"
 * #include "CorrelationProvider.h"
 * #include "QCDProvider.h"
 * #include <iostream>
 *
 * int main() {
 *     ParamId pid{ParameterType::SM, "MASS", 24}; // W boson
 *     
 *     ParameterProvider param_provider(ParameterType::SM);
 *     double w_mass = param_provider(pid);
 *     std::cout << "W mass = " << w_mass << " GeV" << std::endl;
 *     
 *     CorrelationProvider corr_provider;
 *     double correlation = corr_provider(pid, pid, CorrelationProvider::CorrelationType::COMBINED);
 *     std::cout << "W mass correlation with itself = " << correlation << std::endl;
 *     
 *     QCDProvider qcd_provider;
 *     double alpha_s = qcd_provider(AlphasConfig(91.1876, MassType::MSBAR, MassType::MSBAR));
 *     std::cout << "Alpha_s at MZ = " << alpha_s << std::endl;
 *
 *     return 0;
 * }
 * @endcode
 */

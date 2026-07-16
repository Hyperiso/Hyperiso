/**
 * @example examples/freezer_example.cpp
 * @brief Example usage of Freezer to freeze and unfreeze parameters and blocks.
 *
 * @code
 * #include "Freezer.h"
 * #include "Parameters.h"
 * #include <iostream>
 *
 * int main() {
 *     // Assume Parameters are initialized
 *     
 *     ParamId mass_w{ParameterType::SM, "MASS", 24}; // W boson
 *
 *     // Freeze a single parameter
 *     Freezer::freeze(mass_w);
 *     std::cout << "W mass is now frozen." << std::endl;
 *
 *     // Unfreeze a single parameter
 *     Freezer::unfreeze(mass_w);
 *     std::cout << "W mass is now unfrozen." << std::endl;
 *
 *     // Freeze a full block
 *     Freezer::freeze(ParameterType::SM, "GAUGE");
 *     std::cout << "GAUGE block is now frozen." << std::endl;
 *
 *     // Unfreeze the block
 *     Freezer::unfreeze(ParameterType::SM, "GAUGE");
 *     std::cout << "GAUGE block is now unfrozen." << std::endl;
 *
 *     return 0;
 * }
 * @endcode
 */

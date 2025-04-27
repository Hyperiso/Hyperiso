/**
 * @example examples/monitoring_system_example.cpp
 * @brief Example of using monitoring and access classes like HyperisoMaster, MartyAdapter, APIAdapter.
 *
 * @code
 * #include "HyperisoMaster.h"
 * #include "MartyAdapter.h"
 * #include "APIAdapter.h"
 * #include <iostream>
 *
 * int main() {
 *     HyperisoMaster hyperiso;
 *     hyperiso.init("input.lha");
 *
 *     if (hyperiso.check_flag(ExternalFlag::USE_MARTY)) {
 *         MartyAdapter marty;
 *         std::cout << "MARTY Model: " << marty.get_marty_model_name() << std::endl;
 *     }
 *
 *     APIAdapter api;
 *     for (const auto& block : api.get_all_blocks()) {
 *         std::cout << "Available block: " << block << std::endl;
 *     }
 *
 *     return 0;
 * }
 * @endcode
 */

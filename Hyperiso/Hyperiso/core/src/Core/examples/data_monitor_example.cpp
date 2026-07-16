/**
 * @example data_monitor_example.cpp
 * @brief Example usage of APIAdapter to inspect all parameter blocks.
 *
 * @code
 * #include "APIAdapter.h"
 * #include <iostream>
 *
 * int main() {
 *     APIAdapter api;
 *
 *     auto blocks = api.get_all_blocks();
 *     std::cout << "Blocks available:" << std::endl;
 *     for (const auto& block : blocks) {
 *         std::cout << " - " << block << std::endl;
 *     }
 *
 *     return 0;
 * }
 * @endcode
 */

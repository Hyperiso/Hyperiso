/**
 * @example examples/path_provider_example.cpp
 * @brief Example usage of classes implementing IPathProvider (APIAdapter, MartyAdapter).
 *
 * @code
 * #include "APIAdapter.h"
 * #include "MartyAdapter.h"
 * #include <iostream>
 *
 * int main() {
 *     APIAdapter api;
 *     std::cout << "LHA file path: " << api.get_path(APIPath::LHA_PATH) << std::endl;
 *
 *     MartyAdapter marty;
 *     std::cout << "Marty template directory: " << marty.get_path(MartyPath::TEMPLATE_DIR) << std::endl;
 *
 *     return 0;
 * }
 * @endcode
 */

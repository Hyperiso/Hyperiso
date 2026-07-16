#include <iostream>

#include "APIAdapter.h"
#include "BlockProvider.h"
#include "HyperisoMaster.h"
#include "Include.h"
#include "Logger.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    // Dev example: use APIAdapter as a small diagnostic interface for CLI, GUI or notebook tooling.
    HyperisoConfig config;
    HyperisoMaster hyp;
    hyp.init("lha/si_input.flha", config);

    APIAdapter api;

    // Paths are useful to debug which files Hyperiso has loaded.
    std::cout << "LHA path: " << api.get_path(APIPath::LHA_PATH) << "\n";
    std::cout << "Assets root: " << api.get_path(APIPath::ASSETS_ROOT) << "\n";
    std::cout << "User SM params path: " << api.get_path(APIPath::USER_SM_PARAMS) << "\n";

    // Flags are useful to know how the current run was configured.
    std::cout << "HAS_WILSON_INPUT: " << api.check_flag(ExternalFlag::HAS_WILSON_INPUT) << "\n";
    std::cout << "HYP_AS_SM_MARTY: " << api.check_flag(ExternalFlag::HYP_AS_SM_MARTY) << "\n";

    // This can help when a block name is ambiguous, like MASS which can exist in several sections.
    std::cout << "\nPossible types for MASS:\n";
    for (auto type : api.get_type_of_block("MASS")) {
        std::cout << "  " << ParameterTypeMapper::str(type) << "\n";
    }

    // List the first SM blocks and a few values from MASS.
    std::cout << "\nFirst SM blocks:\n";
    std::size_t n = 0;
    for (const auto& block : api.get_blocks_list(ParameterType::SM)) {
        if (n++ >= 10) break;
        std::cout << "  " << block << "\n";
    }

    std::cout << "\nFirst SM::MASS values:\n";
    n = 0;
    for (const auto& [id, value] : api.get_block_infos("MASS", ParameterType::SM)) {
        if (n++ >= 10) break;
        std::cout << "  MASS[" << id << "] = " << value << "\n";
    }

    return 0;
}

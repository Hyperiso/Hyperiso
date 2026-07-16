#include <iostream>

#include "BlockProvider.h"
#include "HyperisoMaster.h"
#include "Include.h"
#include "Logger.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    // All informations for the creation of the HyperisoMaster are provided in the base_core_example.cpp example.
    HyperisoConfig config;
    HyperisoMaster hyp;
    hyp.init("lha/si_input.flha", config);

    // This API allows to search within Hyperiso to get informations on the differents parameters.
    // Can be very helpful when facing a strange result and need to know which parameters are used inside the pipeline.
    BlockProvider block_logger;

    // The parameters are divided in multiple sections:
    //   - ParameterType::SM: Standard Model inputs, such as GAUGE, MASS, SMINPUTS, etc.
    //   - ParameterType::BSM: All model inputs which are not SM. Custom blocks go in this category.
    //   - ParameterType::DECAY: Decays's inputs, such as form factors, etc.
    //   - ParameterType::FLAVOR: Flavor's inputs, such as FCONST, FLIFE, etc.
    //   - ParameterType::WILSON: Scales inputs and Wilson coefficients.
    //   - ParameterType::OBSERVABLE: Experimental observables values.

    std::cout << "\nSMINPUTS:\n";
    // One can get the content of one block using:
    for (const auto& [id, value] : block_logger.get_block(ParameterType::SM, "SMINPUTS")) {
        std::cout << "  " << id << " = " << value << "\n";
    }

    std::cout << "\nWILSON BLOCKS:\n";
    // One can also get all the blocks available within a section.
    for (const auto& block : block_logger.get_all_blocks(ParameterType::WILSON)) {
        std::cout << "  " << block << "\n";
    }

    return 0;
}

#include "HyperisoMaster.h"
#include "ParameterProvider.h"
#include "Logger.h"
#include "config.hpp"
#include <iostream>
#include <cassert>

int main() {

    Logger::getInstance()->setLevel(Logger::LogLevel::DEBUG);

    std::string root_data_file = project_assets_root.data();
    auto mm = HyperisoMaster();

    Config config;
    config.model = Model::SM;
    mm.init(root_data_file + "Test/InputFiles/testInput.flha", config);
    auto sm_provider = ParameterProvider(ParameterType::SM);

    std::cout << "M_Z : " << sm_provider("SMINPUTS", 4) << std::endl;
    assert(std::abs(sm_provider("SMINPUTS", 4) - 9.11876000e+01) < 1e-5);

}
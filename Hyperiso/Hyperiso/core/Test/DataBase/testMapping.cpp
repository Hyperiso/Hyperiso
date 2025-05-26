#include "MappingDataBase.h"
#include "config.hpp"
#include <iostream>
#include <cassert>

int main() {
    std::string root_path = project_root.data();
    std::string jsonFilePath = root_path + "Test/InputFiles/test_mapping.json";
    std::cout << jsonFilePath << std::endl;

    auto dbInstance = MappingDatabase::getInstance("TestInstance", jsonFilePath);
    assert(dbInstance != nullptr);

    const auto& params = dbInstance->getParams();

    assert(params.size() == 2);

    auto param1 = params.find("param1");
    assert(param1 != params.end());
    assert(param1->second.block == "BlockA");
    assert(param1->second.pdgCode == LhaID(1234));

    auto param2 = params.find("param2");
    assert(param2 != params.end());
    assert(param2->second.block == "BlockB");
    assert(param2->second.pdgCode == LhaID(5678));

    auto dbInstanceDuplicate = MappingDatabase::getInstance("TestInstance");
    assert(dbInstance == dbInstanceDuplicate); 

    auto invalidInstance = MappingDatabase::getInstance("InvalidInstance", "invalid_path.json");
    // assert(invalidInstance == nullptr);

    std::cout << "All test succeeded." << std::endl;

    return 0;
}
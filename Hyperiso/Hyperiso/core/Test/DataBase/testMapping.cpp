#include "MappingDatabase.h"
#include "config.hpp"
#include <iostream>
#include <cassert>
#include "JsonParamMappingAdapter.h"
#include "JsonParser.h"

int main() {
    std::string root_path = project_root.data();
    std::string jsonFilePath = root_path + "Test/InputFiles/test_mapping.json";
    std::cout << jsonFilePath << std::endl;
    std::shared_ptr<IParamMappingSource> ips = std::make_shared<JsonParamMappingAdapter>(std::make_shared<JSONParser>());
    std::string name = "TestInstance";
    auto dbInstance = MappingDatabase(name, jsonFilePath, ips);
    // assert(dbInstance != nullptr);

    const auto& params = dbInstance.getParams();

    assert(params.size() == 2);

    auto param1 = params.find("param1");
    assert(param1 != params.end());
    assert(param1->second.block == "BlockA");
    assert(param1->second.code == LhaID(1234));

    auto param2 = params.find("param2");
    assert(param2 != params.end());
    assert(param2->second.block == "BlockB");
    assert(param2->second.code == LhaID(5678));


    auto invalidInstance = MappingDatabase("InvalidInstance", "invalid_path.json", ips);
    // assert(invalidInstance == nullptr);

    std::cout << "All test succeeded." << std::endl;

    return 0;
}
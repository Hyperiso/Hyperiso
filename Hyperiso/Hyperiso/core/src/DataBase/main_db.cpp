#include "MappingDataBase.h"
#include <iostream>

int main() {

    auto map = MappingDatabase::getInstance("SM", "../Models/SM/sm.json");
    const auto& smParams = map->getParams();
    for (const auto& elem : smParams) {
        std::cout << elem.first << std::endl;
        std::cout << elem.second.block << " " << elem.second.pdgCode << std::endl;
    }
    return 0;
}
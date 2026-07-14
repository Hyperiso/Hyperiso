#include "DBManager.h"
#include <iostream>

int main(int argc, char** argv) {

    auto dbm = DBManager();

    const std::string input = argc > 1 ? argv[1] : "Assets/lha/camilia.flha";
    auto node = dbm.read_from_file(input);

    node->printJSON();
    return 0;
}
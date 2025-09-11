#include "DBManager.h"
#include <iostream>

int main() {

    auto dbm = DBManager();

    auto node = dbm.read_from_file("/home/theo/hyperiso/Assets/lha/camilia.flha");

    node->printJSON();
    return 0;
}
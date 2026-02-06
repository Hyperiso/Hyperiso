#include "DBManager.h"
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>

namespace fs = std::filesystem;

int main() {
    std::cout << "== Running UNIT tests for DBManager ==\n";

    {
        DBManager db;
        bool threw = false;
        try {
            (void)db.read_from_file("does_not_exist.json");
        } catch (...) {
            threw = true;
        }
        assert(threw);
    }

    {
        DBManager db;
        auto node = std::make_shared<DBNode>();
        node->set("whatever", "somekey");

        bool threw = false;
        try {
            db.write_to_file("tmp_unit.json", node);
        } catch (...) {
            threw = true;
        }
        assert(threw);


        node->set(1, "SMINPUTS");
        db.write_to_file("tmp_unit.json", node);
        assert(fs::exists("tmp_unit.json"));
        fs::remove("tmp_unit.json");
    }

    {
        DBManager db;
        db.add_lha_prototype("MYBLOCK", 2, 1);
        db.add_lha_prototype("MYBLOCK", 2, 1); 
    }

    std::cout << "\n All DBManager unit tests passed!\n";
    return 0;
}

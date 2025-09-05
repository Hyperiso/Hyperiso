// test_DBManager_unit.cpp
#include "DBManager.h"
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>

namespace fs = std::filesystem;

int main() {
    std::cout << "== Running UNIT tests for DBManager ==\n";

    // 1) read_from_file sur fichier manquant -> doit throw
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

    // 2) write_to_file exige "SMINPUTS"
    {
        DBManager db;
        auto node = std::make_shared<Node>();
        node->set("whatever", "somekey");

        bool threw = false;
        try {
            db.write_to_file("tmp_unit.json", node);
        } catch (...) {
            threw = true;
        }
        assert(threw); // manque SMINPUTS -> OK

        // Ajoute la clé requise et ré-essaie
        node->set(1, "SMINPUTS");
        db.write_to_file("tmp_unit.json", node);
        assert(fs::exists("tmp_unit.json"));
        fs::remove("tmp_unit.json");
    }

    // 3) add_lha_prototype idempotent (même prototype = warning mais pas d'exception)
    {
        DBManager db;
        db.add_lha_prototype("MYBLOCK", 2, 1);
        db.add_lha_prototype("MYBLOCK", 2, 1); // ne doit pas throw
    }

    std::cout << "\n✅ All DBManager unit tests passed!\n";
    return 0;
}

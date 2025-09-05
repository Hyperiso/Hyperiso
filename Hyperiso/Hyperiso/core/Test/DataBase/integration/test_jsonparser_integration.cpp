// test_jsonparser_integration.cpp
#include "JsonParser.h"
#include "DBNode.h"
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cmath>

namespace fs = std::filesystem;

static void create_json_file(const fs::path& path) {
    std::ofstream f(path);
    f << R"({"model":{"type":"SM","MZ":91.1876},"parameters":{"alpha_s":0.118,"top_mass":173.1}})";
}

int main() {
    std::cout << "== Running INTEGRATION tests for JSONParser ==\n";

    const fs::path in = "jp_in.json";
    const fs::path out = "jp_out.json";
    create_json_file(in);

    JSONParser p;

    // lecture
    auto root = p.readFromFile(in);
    assert(std::get<BlockName>(root->get("model","type")) == "SM");
    assert(std::abs(std::get<double>(root->get("parameters","alpha_s")) - 0.118) < 1e-12);

    // modifs
    root->set("2HDM", "model", "type");
    root->set(125.0, "model", "MH");

    // écriture
    p.writeToFile(out.string(), root);

    // relecture
    auto back = p.readFromFile(out.string());
    assert(std::get<BlockName>(back->get("model","type")) == "2HDM");
    assert(std::abs(std::get<double>(back->get("model","MH")) - 125.0) < 1e-12);

    fs::remove(in);
    fs::remove(out);

    std::cout << "\n✅ JSONParser integration tests passed!\n";
    return 0;
}

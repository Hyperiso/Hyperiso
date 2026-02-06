// test_DBManager_integration.cpp
#include "DBManager.h"
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cmath>

namespace fs = std::filesystem;

static void create_test_json_file(const fs::path& path) {
    std::ofstream file(path);
    file << R"({
        "model": { "type": "SM", "MZ": 91.1876 },
        "parameters": { "alpha_s": 0.118, "top_mass": 173.1 }
    })";
}

static void create_test_yaml_file(const fs::path& path) {
    std::ofstream file(path);
    file << R"(model:
  type: SM
  MZ: 91.1876
parameters:
  alpha_s: 0.118
  top_mass: 173.1
)";
}

static void create_test_flha_file(const fs::path& path) {
    std::ofstream file(path);
    file << R"(Block MASS  # Mass spectrum
   25  1.25000000E+02  # h0
   35  2.00000000E+02  # H0
Block GAUGE
   1   3.57300000E-01  # g1
   2   6.46400000E-01  # g2
   3   1.21700000E+00  # g3
)";
}

static void test_read_modify_write_and_reload_json() {
    std::cout << "\n-- test_read_modify_write_and_reload_json --" << std::endl;

    const fs::path in = "test_in.json";
    const fs::path out = "test_out.json";
    create_test_json_file(in);

    DBManager db;
    auto node = db.read_from_file(in);

    assert(std::get<BlockName>(node->get("model", "type")) == "SM");
    assert(std::abs(std::get<double>(node->get("parameters", "alpha_s")) - 0.118) < 1e-6);

    node->set("2HDM", "model", "type");
    node->set(125.0, "model", "MH");

    node->set(1, "SMINPUTS");

    db.write_to_file(out, node);

    auto reloaded = db.read_from_file(out);
    assert(std::get<BlockName>(reloaded->get("model", "type")) == "2HDM");
    assert(std::abs(std::get<double>(reloaded->get("model", "MH")) - 125.0) < 1e-6);

    fs::remove(in);
    fs::remove(out);
}

static void test_yaml_read_modify_write_reload() {
    std::cout << "\n-- test_yaml_read_modify_write_reload --" << std::endl;

    const fs::path input_file = "test_in.yaml";
    const fs::path output_file = "test_out.yaml";

    create_test_yaml_file(input_file);

    DBManager db;
    auto node = db.read_from_file(input_file);

    assert(std::get<BlockName>(node->get("model", "type")) == "SM");
    assert(std::abs(std::get<double>(node->get("parameters", "top_mass")) - 173.1) < 1e-6);

    node->set("YAML_MODE", "meta", "type");
    node->set(999, "meta", "version");

    // clé requise par DBManager::write_to_file
    node->set(1, "SMINPUTS");

    // plus de redirection cout: le YAML est écrit dans output_file
    db.write_to_file(output_file, node);

    auto reloaded = db.read_from_file(output_file);
    assert(std::get<BlockName>(reloaded->get("meta", "type")) == "YAML_MODE");
    assert(std::get<int>(reloaded->get("meta", "version")) == 999);

    fs::remove(input_file);
    fs::remove(output_file);
}

static double as_double_or_string(const DBNode::Value& v) {
    if (std::holds_alternative<double>(v)) {
        return std::get<double>(v);
    }
    if (std::holds_alternative<BlockName>(v)) {
        return std::stod(std::get<BlockName>(v));
    }
    // Si jamais c'est un noeud, on lève
    throw std::runtime_error("Expected numeric (double/string), got non-scalar");
}

static double as_number(const DBNode::Value& v) {
    if (std::holds_alternative<double>(v)) {
        return std::get<double>(v);
    }
    if (std::holds_alternative<BlockName>(v)) {
        return std::stod(std::get<BlockName>(v));
    }
    if (std::holds_alternative<std::shared_ptr<DBNode>>(v)) {
        auto sub = std::get<std::shared_ptr<DBNode>>(v);
        // lecture directe de central_value (double ou string)
        auto cv = sub->get("central_value");
        if (std::holds_alternative<double>(cv)) {
            return std::get<double>(cv);
        }
        if (std::holds_alternative<BlockName>(cv)) {
            return std::stod(std::get<BlockName>(cv));
        }
    }
    throw std::runtime_error("Expected numeric (double/string), got non-scalar");
}

static void test_flha_read_and_extract_mass() {
    std::cout << "\n-- test_flha_read_and_extract_mass --" << std::endl;

    const fs::path path = "test_model.flha";
    create_test_flha_file(path);

    DBManager db;
    // Ton fichier a 2 colonnes => garde ces prototypes
    db.add_lha_prototype("MASS",  2, 1);
    db.add_lha_prototype("GAUGE", 2, 1);

    auto node = db.read_from_file(path);

    assert(node->contains("MASS"));
    assert(node->contains("GAUGE"));

    auto gauge_group = node->getGroup({"GAUGE"});
    assert(std::abs(as_number(gauge_group["1"]) - 0.3573) < 1e-4);
    assert(std::abs(as_number(gauge_group["3"]) - 1.217)  < 1e-4);

    auto mass_group = node->getGroup({"MASS"});
    assert(std::abs(as_number(mass_group["25"]) - 125.0) < 1e-2);
    assert(std::abs(as_number(mass_group["35"]) - 200.0) < 1e-2);

    fs::remove(path);
}

int main() {
    std::cout << "== Running INTEGRATION test for DBManager ==\n";
    test_read_modify_write_and_reload_json();
    test_yaml_read_modify_write_reload();
    test_flha_read_and_extract_mass();
    std::cout << "\n DBManager integration test passed!\n";
    return 0;
}

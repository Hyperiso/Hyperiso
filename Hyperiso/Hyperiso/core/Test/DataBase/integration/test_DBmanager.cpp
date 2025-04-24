#include "DBManager.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <filesystem>

namespace fs = std::filesystem;

void create_test_json_file(const fs::path& path) {
    std::ofstream file(path);
    file << R"({
        "model": {
            "type": "SM",
            "MZ": 91.1876
        },
        "parameters": {
            "alpha_s": 0.118,
            "top_mass": 173.1
        }
    })";
    file.close();
}

void test_read_modify_write_and_reload() {
    std::cout << "\n-- test_read_modify_write_and_reload --" << std::endl;

    fs::path input_file = "test_input.json";
    fs::path output_file = "test_output.json";

    create_test_json_file(input_file);

    DBManager db;
    auto node = db.read_from_file(input_file);

    // Vérif lecture
    assert(std::get<std::string>(node->get("model", "type")) == "SM");
    assert(std::abs(std::get<double>(node->get("parameters", "alpha_s")) - 0.118) < 1e-6);

    // Modification
    node->set("2HDM", "model", "type");
    node->set(125.0, "model", "MH");

    db.write_to_file(output_file, node);

    // Relecture
    auto reloaded = db.read_from_file(output_file);
    assert(std::get<std::string>(reloaded->get("model", "type")) == "2HDM");
    assert(std::abs(std::get<double>(reloaded->get("model", "MH")) - 125.0) < 1e-6);

    // Nettoyage fichiers temporaires
    fs::remove(input_file);
    fs::remove(output_file);
}

void create_test_yaml_file(const fs::path& path) {
    std::ofstream file(path);
    file << R"(
model:
  type: SM
  MZ: 91.1876
parameters:
  alpha_s: 0.118
  top_mass: 173.1
)";
    file.close();
}

void test_yaml_read_modify_write_reload() {
    std::cout << "\n-- test_yaml_read_modify_write_reload --" << std::endl;

    fs::path input_file = "test_input.yaml";
    fs::path output_file = "test_output.yaml";

    create_test_yaml_file(input_file);

    DBManager db;
    auto node = db.read_from_file(input_file);

    // Vérif lecture
    assert(std::get<std::string>(node->get("model", "type")) == "SM");
    assert(std::abs(std::get<double>(node->get("parameters", "top_mass")) - 173.1) < 1e-6);

    // Modifications
    node->set("YAML_MODE", "meta", "type");
    node->set(999, "meta", "version");

    db.write_to_file(output_file, node);

    auto reloaded = db.read_from_file(output_file);
    assert(std::get<std::string>(reloaded->get("meta", "type")) == "YAML_MODE");
    assert(std::get<int>(reloaded->get("meta", "version")) == 999);

    fs::remove(input_file);
    fs::remove(output_file);
}

void create_test_flha_file(const fs::path& path) {
    std::ofstream file(path);
    file << R"(
Block MASS  # Mass spectrum
   25  1.25000000E+02  # h0
   35  2.00000000E+02  # H0
Block GAUGE
   1   3.57300000E-01  # g1
   2   6.46400000E-01  # g2
   3   1.21700000E+00  # g3
)";
    file.close();
}

void test_flha_read_and_extract_mass() {
    std::cout << "\n-- test_flha_read_and_extract_mass --" << std::endl;

    fs::path flha_path = "test_model.flha";
    create_test_flha_file(flha_path);

    DBManager db;
    db.add_lha_prototype("MASS", 2, 1);
    db.add_lha_prototype("GAUGE", 2, 1);

    auto node = db.read_from_file(flha_path);

    assert(node->contains("MASS"));
    assert(node->contains("GAUGE"));

    auto gauge_group = node->getGroup({"GAUGE"});
    assert(std::abs(std::get<double>(gauge_group["1"]) - 0.3573) < 1e-4);
    assert(std::abs(std::get<double>(gauge_group["3"]) - 1.217) < 1e-4);

    auto mass_group = node->getGroup({"MASS"});
    assert(std::abs(std::get<double>(mass_group["25"]) - 125.0) < 1e-2);
    assert(std::abs(std::get<double>(mass_group["35"]) - 200.0) < 1e-2);

    fs::remove(flha_path);
}

int main() {
    std::cout << "== Running INTEGRATION test for DBManager ==\n";

    test_read_modify_write_and_reload();
    test_yaml_read_modify_write_reload();
    test_flha_read_and_extract_mass();
    
    std::cout << "\n✅ DBManager integration test passed!\n" << std::endl;
    return 0;
}

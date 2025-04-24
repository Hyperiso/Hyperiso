#include "DBNode.h"
#include <iostream>
#include <cassert>
#include <sstream>

void test_complex_hierarchy_and_output() {
    std::cout << "\n-- test_complex_hierarchy_and_output --" << std::endl;

    Node root;

    // Section "model"
    root.set(std::string("SM"), "model", "type");
    root.set(true, "model", "use_QED");
    root.set(91.1876, "model", "MZ");

    // Section "parameters"
    root.set(0.118, "parameters", "alpha_s");
    root.set(173.1, "parameters", "top_mass");

    // Section "metadata"
    root.set(42, "metadata", "run_id");
    root.set("integration test", "metadata", "description");

    // Section "observables" as list of Node
    std::vector<std::shared_ptr<Node>> obs_list;

    auto obs1 = std::make_shared<Node>();
    obs1->set("BR(B→X_s γ)", "");
    auto obs2 = std::make_shared<Node>();
    obs2->set("R(D*)", "");

    obs_list.push_back(obs1);
    obs_list.push_back(obs2);
    root.set(obs_list, "observables");

    // Validate values
    assert(std::get<std::string>(root.get("model", "type")) == "SM");
    assert(std::abs(std::get<double>(root.get("model", "MZ")) - 91.1876) < 1e-6);
    assert(std::get<bool>(root.get("model", "use_QED")) == true);
    assert(std::get<int>(root.get("metadata", "run_id")) == 42);

    // Print output as JSON
    std::cout << "\n--- JSON Output ---\n";
    root.printJSON();
    std::cout << std::endl;

    // Print output as YAML
    std::cout << "\n--- YAML Output ---\n";
    root.printYAML();
    std::cout << std::endl;
}

void test_json_stream_output() {
    std::cout << "\n-- test_json_stream_output --" << std::endl;

    Node root;
    root.set("streaming test", "info");
    root.set(123, "version");

    std::stringstream ss;
    root.printJSONToStream(ss);

    std::string result = ss.str();
    std::cout << result << std::endl;

    assert(result.find("\"info\"") != std::string::npos);
    assert(result.find("streaming test") != std::string::npos);
    assert(result.find("123") != std::string::npos);
}

int main() {
    std::cout << "== Running INTEGRATION tests for Node ==\n";

    test_complex_hierarchy_and_output();
    test_json_stream_output();

    std::cout << "\n✅ All Node integration tests passed!\n" << std::endl;
    return 0;
}

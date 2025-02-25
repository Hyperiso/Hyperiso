#include "Parser.h"
#include "DBNode.h"
#include <memory>
#include <string>
#include <iostream>
#include <fstream>

int main_YAML() {
    YAMLParser yamlParser;

    std::cout << "\n🔹 Parsing YAML from string:\n";
    std::string yamlInput = R"(
name: Bob
age: 25
15: 
  truc: 5
  machin: 1.5
location:
  city: London
  country: UK
skills:
  - JavaScript
  - Rust
  - Go
)";
    auto nodeYaml = yamlParser.parse(yamlInput);
    nodeYaml->printJSON();

    try {
        auto name = std::get<std::string>(nodeYaml->get("name"));
        std::cout << "✅ Name: " << name << "\n";
    } catch (const std::exception& e) {
        std::cerr << "❌ Error: " << e.what() << "\n";
    }

    try {
        auto city = std::get<std::string>(nodeYaml->get("location", "city"));
        std::cout << "✅ City: " << city << "\n";
    } catch (const std::exception& e) {
        std::cerr << "❌ Error: " << e.what() << "\n";
    }

    try {
        auto skill1 = std::get<std::string>(nodeYaml->get("skills", "0"));
        auto skill2 = std::get<std::string>(nodeYaml->get("skills", "1"));
        auto skill3 = std::get<std::string>(nodeYaml->get("skills", "2"));
        std::cout << "✅ Skills: " << skill1 << ", " << skill2 << ", " << skill3 << "\n";
    } catch (const std::exception& e) {
        std::cerr << "❌ Error: " << e.what() << "\n";
    }

    try {
        auto group = nodeYaml->getGroup({"location"});
        std::cout << "✅ Group 'location':\n";
        for (const auto& [key, val] : group) {
            std::cout << "- " << key << " : ";
            if (std::holds_alternative<std::string>(val)) {
                std::cout << std::get<std::string>(val);
            } else if (std::holds_alternative<int>(val)) {
                std::cout << std::get<int>(val);
            }
            std::cout << "\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "❌ Error: " << e.what() << "\n";
    }

    std::map<std::string, Node::Value> newGroup = {
        {"hobby", "cycling"},
        {"job", "engineer"}
    };
    nodeYaml->setGroup({"extra"}, newGroup);
    std::cout << "\n✅ Structure YAML after setGroup:\n";
    nodeYaml->printJSON();

    std::cout << "\n🔹 Writing YAML to file...\n";
    std::ofstream file("test_output.yaml");
    if (file.is_open()) {
        nodeYaml->printYAML();
        file.close();
        std::cout << "✅ YAML File written successfully.\n";
    } else {
        std::cerr << "❌ Error writing YAML file.\n";
    }

    std::cout << "\n🔹 Checking written YAML content:\n";
    std::ifstream fileCheck("test_output.yaml");
    if (fileCheck.is_open()) {
        std::ostringstream ossCheck;
        ossCheck << fileCheck.rdbuf();
        std::cout << ossCheck.str() << "\n";
        fileCheck.close();
    } else {
        std::cerr << "❌ Error checking YAML file content.\n";
    }

    std::cout << "\n🔹 Reading YAML from file...\n";
    std::ifstream fileRead("test_output.yaml");
    if (fileRead.is_open()) {
        std::ostringstream oss;
        oss << fileRead.rdbuf();
        std::cout << "✅ YAML File read before parsing:\n" << oss.str() << "\n";
        auto parsedFileNode = yamlParser.parse(oss.str());
        std::cout << "✅ YAML File read successfully:\n";
        parsedFileNode->printJSON();
        fileRead.close();
    } else {
        std::cerr << "❌ Error reading YAML file.\n";
    }

    return 0;
}

int main() {
    JSONParser jsonParser;
    YAMLParser yamlParser;

    std::string jsonInput = R"({
        "name": "Alice",
        "age": 1.2e-1,
        "details": {
            "city": "Paris",
            "skills": ["C++", "Python", "YAML"]
        }
    })";
    std::cout << "🔹 Parsing JSON:\n";
    auto nodeJson = jsonParser.parse(jsonInput);
    nodeJson->printJSON();

    std::string yamlInput = R"(
name: Alice
age: 30
details:
  city: Paris
  skills:
    - C++
    - Python
    - YAML
)";
    std::cout << "\n🔹 Parsing YAML:\n";
    auto nodeYaml = yamlParser.parse(yamlInput);
    nodeYaml->printJSON();

    std::cout << "\n🔹 Testing Node Manipulations:\n";
    auto root = std::make_shared<Node>();
    root->set("Alice", "name");
    root->set(1.5e-2, "age");
    root->set(true, "isAdmin");
    root->set("Paris", "address", "city");
    root->set(75001, "address", "zip");

    auto listNode = std::make_shared<Node>();
    listNode->set("C++", "0");
    listNode->set("Python", "1");
    listNode->set("YAML", "2");
    root->set(listNode, "skills");

    std::cout << "✅ Generated JSON from Node:\n";
    root->printJSON();

    try {
        auto name = std::get<std::string>(root->get("name"));
        std::cout << "✅ Name: " << name << "\n";
    } catch (const std::exception& e) {
        std::cerr << "❌ Error: " << e.what() << "\n";
    }

    try {
        auto city = std::get<std::string>(root->get("address", "city"));
        std::cout << "✅ Address city: " << city << "\n";
    } catch (const std::exception& e) {
        std::cerr << "❌ Error: " << e.what() << "\n";
    }

    try {
        auto skill1 = std::get<std::string>(root->get("skills", "0"));
        auto skill2 = std::get<std::string>(root->get("skills", "1"));
        auto skill3 = std::get<std::string>(root->get("skills", "2"));
        std::cout << "✅ Skills: " << skill1 << ", " << skill2 << ", " << skill3 << "\n";
    } catch (const std::exception& e) {
        std::cerr << "❌ Error: " << e.what() << "\n";
    }

    try {
        auto group = root->getGroup({"address"});
        std::cout << "✅ Group 'address':\n";
        for (const auto& [key, val] : group) {
            std::cout << "- " << key << " : ";
            if (std::holds_alternative<std::string>(val)) {
                std::cout << std::get<std::string>(val);
            } else if (std::holds_alternative<int>(val)) {
                std::cout << std::get<int>(val);
            }
            std::cout << "\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "❌ Error: " << e.what() << "\n";
    }

    std::map<std::string, Node::Value> newGroup = {
        {"newKey1", "newValue1"},
        {"newKey2", 12345},
        {"newKey3", false},
    };
    root->setGroup({"newSection"}, newGroup);
    std::cout << "\n✅ Structure JSON after setGroup:\n";
    root->printJSON();

    std::cout << "\n🔹 Writing to file...\n";
    std::ofstream file("test_output.json");
    if (file.is_open()) {
        root->printJSONToStream(file);
        file.close();
        std::cout << "✅ File written successfully.\n";
    } else {
        std::cerr << "❌ Error writing file.\n";
    }

    std::cout << "\n🔹 Reading from file...\n";
    std::ifstream fileRead("test_output.json");
    if (fileRead.is_open()) {
        std::ostringstream oss;
        oss << fileRead.rdbuf();
        auto parsedFileNode = jsonParser.parse(oss.str());
        std::cout << "✅ File read successfully:\n";
        parsedFileNode->printJSON();
        fileRead.close();
    } else {
        std::cerr << "❌ Error reading file.\n";
    }




    main_YAML();
    return 0;
}

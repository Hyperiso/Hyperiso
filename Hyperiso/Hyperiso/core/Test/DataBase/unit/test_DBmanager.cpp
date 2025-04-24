#include "DBManager.h"
#include <iostream>
#include <cassert>

// === Mock parsers ===

class MockParser : public IParser {
public:
    std::shared_ptr<Node> readFromFile(const fs::path& path) override {
        auto node = std::make_shared<Node>();
        node->set("mock_value", "key");
        return node;
    }

    void writeToFile(const fs::path& path, std::shared_ptr<Node> root) override {
        std::cout << "[MockParser] Pretending to write to " << path << std::endl;
        assert(root->contains("key"));  // we expect the mock node to have "key"
    }
};

// Override ParserFactory temporarily
namespace {
    struct ScopedParserOverride {
        ScopedParserOverride() {
            ParserFactory::createParser = [](ParserFactory::Type type) {
                return std::make_shared<MockParser>();
            };
        }
    };
}

// === Unit tests ===

void test_deduce_parser_type() {
    std::cout << "\n-- test_deduce_parser_type --" << std::endl;
    DBManager db;
    
    assert(db.read_from_file("test.json")->get_keys().size() > 0);
    assert(db.read_from_file("file.flha")->get("key") == "mock_value");
}

void test_add_lha_prototype_duplicate() {
    std::cout << "\n-- test_add_lha_prototype_duplicate --" << std::endl;
    DBManager db;
    db.add_lha_prototype("MYBLOCK", 2);

    // Répétition volontaire du même prototype
    db.add_lha_prototype("MYBLOCK", 2);  // Should log a warning, not throw

    // Prototypage différent : devrait lever une erreur
    try {
        db.add_lha_prototype("MYBLOCK", 3);
        assert(false && "Expected exception for conflicting prototype");
    } catch (const std::exception& e) {
        std::cout << "Caught expected exception: " << e.what() << std::endl;
    }
}

int main() {
    std::cout << "== Running UNIT tests for DBManager ==\n";

    // Override parser factory temporarily
    ScopedParserOverride _;

    test_deduce_parser_type();
    test_add_lha_prototype_duplicate();

    std::cout << "\n✅ All DBManager unit tests passed!\n" << std::endl;
    return 0;
}

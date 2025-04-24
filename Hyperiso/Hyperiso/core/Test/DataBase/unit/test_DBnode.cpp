#include "DBNode.h"
#include <iostream>
#include <cassert>
#include <stdexcept>

void test_simple_set_and_get() {
    std::cout << "\n-- test_simple_set_and_get --" << std::endl;

    Node node;
    node.set(42, "intVal");
    node.set(3.14, "doubleVal");
    node.set(std::string("hello"), "strVal");
    node.set(true, "boolVal");

    assert(std::get<int>(node.get("intVal")) == 42);
    assert(std::abs(std::get<double>(node.get("doubleVal")) - 3.14) < 1e-6);
    assert(std::get<std::string>(node.get("strVal")) == "hello");
    assert(std::get<bool>(node.get("boolVal")) == true);
}

void test_nested_set_and_get() {
    std::cout << "\n-- test_nested_set_and_get --" << std::endl;

    Node node;
    node.set(123, "level1", "level2", "value");

    assert(std::get<int>(node.get("level1", "level2", "value")) == 123);
}

void test_get_keys_and_contains() {
    std::cout << "\n-- test_get_keys_and_contains --" << std::endl;

    Node node;
    node.set("data", "key1");
    node.set(100, "key2");

    auto keys = node.get_keys();
    assert(keys.size() == 2);
    assert(node.contains("key1"));
    assert(node.contains("key2"));
    assert(!node.contains("missing"));
}

void test_setGroup_and_getGroup() {
    std::cout << "\n-- test_setGroup_and_getGroup --" << std::endl;

    Node node;
    std::map<std::string, Node::Value> group = {
        {"a", 1},
        {"b", std::string("x")},
        {"c", 2.2}
    };

    node.setGroup({"group1", "sub"}, group);

    auto g = node.getGroup({"group1", "sub"});
    assert(std::get<int>(g["a"]) == 1);
    assert(std::get<std::string>(g["b"]) == "x");
    assert(std::abs(std::get<double>(g["c"]) - 2.2) < 1e-6);
}

void test_get_invalid_path_throws() {
    std::cout << "\n-- test_get_invalid_path_throws --" << std::endl;

    Node node;
    node.set(42, "valid");

    try {
        node.get("nonexistent");
        assert(false && "Expected exception not thrown");
    } catch (const std::runtime_error& e) {
        std::cout << "Caught expected: " << e.what() << std::endl;
    }
}

int main() {
    std::cout << "== Running UNIT tests for Node ==\n";

    test_simple_set_and_get();
    test_nested_set_and_get();
    test_get_keys_and_contains();
    test_setGroup_and_getGroup();
    test_get_invalid_path_throws();

    std::cout << "\n✅ All Node unit tests passed!\n" << std::endl;
    return 0;
}

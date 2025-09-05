// test_jsonparser_unit.cpp
#include "JsonParser.h"
#include "DBNode.h"
#include <cassert>
#include <iostream>
#include <sstream>
#include <cmath>

static bool contains_all(const std::string& s, std::initializer_list<const char*> needles) {
    for (auto* n : needles) if (s.find(n) == std::string::npos) return false;
    return true;
}

int main() {
    std::cout << "== Running UNIT tests for JSONParser ==\n";
    JSONParser p;

    // 1) Objet simple + imbriqué + primitifs
    {
        const std::string json =
            R"({"model":{"type":"SM","MZ":91.1876},"flags":{"use_QED":true},"n":-12,"x":1.23e-2})";
        auto root = p.parse(json);
        assert(std::get<BlockName>(root->get("model","type")) == "SM");
        assert(std::abs(std::get<double>(root->get("model","MZ")) - 91.1876) < 1e-9);
        assert(std::get<bool>(root->get("flags","use_QED")) == true);
        assert(std::abs(std::get<double>(root->get("n")) + 12.0) < 1e-12); // nombres => double
        assert(std::abs(std::get<double>(root->get("x")) - 1.23e-2) < 1e-12);
    }

    // 2) Tableaux => parsés comme sous-noeuds avec clés "0","1",...
    {
        const std::string json = R"({"arr":[1,2,3],"strs":["a","b"],"objs":[{"a":1},{"b":2}]})";
        auto root = p.parse(json);

        auto arr = std::get<std::shared_ptr<Node>>(root->get("arr"));
        assert(std::abs(std::get<double>(arr->get("0")) - 1.0) < 1e-9);
        assert(std::abs(std::get<double>(arr->get("1")) - 2.0) < 1e-9);
        assert(std::abs(std::get<double>(arr->get("2")) - 3.0) < 1e-9);

        auto strs = std::get<std::shared_ptr<Node>>(root->get("strs"));
        assert(std::get<BlockName>(strs->get("0")) == "a");
        assert(std::get<BlockName>(strs->get("1")) == "b");

        auto objs = std::get<std::shared_ptr<Node>>(root->get("objs"));
        auto o0 = std::get<std::shared_ptr<Node>>(objs->get("0"));
        auto o1 = std::get<std::shared_ptr<Node>>(objs->get("1"));
        assert(std::abs(std::get<double>(o0->get("a")) - 1.0) < 1e-9);
        assert(std::abs(std::get<double>(o1->get("b")) - 2.0) < 1e-9);
    }

    // 3) Erreurs : JSON invalide
    {
        bool threw = false;
        try {
            (void)p.parse(R"({model: "SM"})"); // clé non quotée
        } catch (...) { threw = true; }
        assert(threw);

        threw = false;
        try {
            (void)p.parse(R"({"x": tru})"); // bool tronqué
        } catch (...) { threw = true; }
        assert(threw);

        threw = false;
        try {
            // nombre invalide (double 'e')
            (void)p.parse(R"({"x": 1.2ee5})");
        } catch (...) { threw = true; }
        assert(threw);
    }

    // 4) printJSONToStream sur arbre parsé
    {
        const std::string json = R"({"a":1,"b":"x"})";
        auto root = p.parse(json);
        std::stringstream ss;
        root->printJSONToStream(ss);
        auto out = ss.str();
        assert(contains_all(out, {"\"a\"", "1", "\"b\"", "\"x\""}));
    }

    std::cout << "\n✅ All JSONParser unit tests passed!\n";
    return 0;
}

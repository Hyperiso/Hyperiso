// test_lha_blocks_unit.cpp
#include <cassert>
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include <algorithm>
#include <sstream>

#include "lha_blocks.h"
#include "LhaBlockPrototype.h"
#include "lha_elements.h"
#include "DBNode.h"

static bool contains_all(const std::string& s, const std::vector<std::string>& needles) {
    for (auto& n : needles) if (s.find(n) == std::string::npos) return false;
    return true;
}

static double as_number_from_value(const Node::Value& v) {
    if (std::holds_alternative<double>(v)) return std::get<double>(v);
    if (std::holds_alternative<BlockName>(v)) return std::stod(std::get<BlockName>(v));
    if (std::holds_alternative<std::shared_ptr<Node>>(v)) {
        auto sub = std::get<std::shared_ptr<Node>>(v);
        auto cv = sub->get("central_value");
        if (std::holds_alternative<double>(cv)) return std::get<double>(cv);
        if (std::holds_alternative<BlockName>(cv)) return std::stod(std::get<BlockName>(cv));
    }
    throw std::runtime_error("not a numeric-compatible value");
}

int main() {
    std::cout << "== Running UNIT tests for LhaBlock ==\n";

    // 1) Bloc sans échelle globale : SMINPUTS
    {
        LhaBlock blk(SMINPUTS); // itemCount=2, valueIdx=1
        blk.addElement({"1", "0.118"});
        blk.addElement({"2", "173.1"});

        // entries
        auto entries = blk.getEntries();
        assert(entries && entries->size() == 2);

        // hasElement / get
        assert(blk.hasElement(LhaID({1})));
        assert(blk.hasElement(LhaID({2})));
        assert(!blk.hasElement(LhaID({3})));

        auto e1 = blk.get(LhaID({1}));
        auto e2 = blk.get(LhaID({2}));
        assert(e1 && e2);

        // prototype renvoyé
        auto p = blk.getPrototype();
        assert(p.blockName == "SMINPUTS");

        // toString : contient au moins l’en-tête et les valeurs
        auto s = blk.toString();
        assert(contains_all(s, {"Block SMINPUTS", "0.118", "173.1"}));

        // toDBNode : pas de clé "scale", group SMINPUTS présent
        auto n = blk.toDBNode();
        assert(n && !n->contains("scale"));
        auto g = n->getGroup({"SMINPUTS"});
        // Deux éléments, valeurs récupérables
        assert(g.size() == 2);
        // On ne dépend pas du format de la clé (LhaID.to_string()), on lit juste les valeurs
        int seen = 0;
        for (auto& [k, v] : g) {
            double x = as_number_from_value(v);
            if (std::abs(x - 0.118) < 1e-9) ++seen;
            if (std::abs(x - 173.1) < 1e-9) ++seen;
        }
        assert(seen == 2);
    }

    // 2) Bloc avec échelle globale : HMIX (globalScale=true)
    {
        LhaBlock blk(HMIX); // itemCount=3, valueIdx=2, globalScale=true
        // Q (=col 0) dans chaque ligne
        blk.addElement({"1000.0", "1", "0.5"});
        blk.addElement({"1000.0", "2", "1.2"});

        // toDBNode : "scale" présent et = première échelle
        auto n = blk.toDBNode();
        assert(n && n->contains("scale"));
        assert(std::abs(std::get<double>(n->get("scale")) - 1000.0) < 1e-12);

        // group HMIX a 2 éléments, valeurs ok
        auto g = n->getGroup({"HMIX"});
        assert(g.size() == 2);
        int seen = 0;
        for (auto& [k, v] : g) {
            double x = as_number_from_value(v);
            if (std::abs(x - 0.5) < 1e-12) ++seen;
            if (std::abs(x - 1.2) < 1e-12) ++seen;
        }
        assert(seen == 2);
    }

    // 3) readData() ignore les lignes vides
    {
        LhaBlock blk(SMINPUTS);
        std::vector<std::vector<std::string>> lines = {
            {"1", "0.119"}, {}, {"2", "80.379"}
        };
        blk.readData(lines);
        assert(blk.getEntries()->size() == 2);
    }

    // 4) get() renvoie nullptr si absent
    {
        LhaBlock blk(MASS);
        blk.addElement({"25", "125.0"});
        assert(blk.get(LhaID({24})) == nullptr);
    }

    std::cout << "\n✅ All LhaBlock unit tests passed!\n";
    return 0;
}

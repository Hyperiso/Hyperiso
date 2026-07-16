// test_dbnode_unit.cpp
#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>
#include <memory>
#include <map>
#include <algorithm>
#include <cmath>
#include <type_traits>

#include "DBNode.h"

static bool contains_all(const std::string& s, const std::vector<std::string>& needles) {
    for (auto& n : needles) if (s.find(n) == std::string::npos) return false;
    return true;
}

int main() {
    std::cout << "== Running UNIT tests for DBNode ==\n";

    // 1) set/get de tous les types supportés
    {
        DBNode n;
        n.set(42, "intVal");
        n.set(3.14, "doubleVal");
        n.set(std::string("hello"), "strVal");
        n.set(true, "boolVal");

        // shared_ptr<DBNode>
        auto sub = std::make_shared<DBNode>();
        sub->set("leaf", "k");
        n.set(sub, "subNode");

        std::vector<std::shared_ptr<DBNode>> vec;
        auto a = std::make_shared<DBNode>(); a->set("A", "");
        auto b = std::make_shared<DBNode>(); b->set("B", "");
        vec.push_back(a); vec.push_back(b);
        n.set(vec, "list");

        assert(std::get<int>(n.get("intVal")) == 42);
        assert(std::abs(std::get<double>(n.get("doubleVal")) - 3.14) < 1e-9);
        assert(std::get<BlockName>(n.get("strVal")) == "hello");
        assert(std::get<bool>(n.get("boolVal")) == true);

        auto sp = std::get<std::shared_ptr<DBNode>>(n.get("subNode"));
        assert(sp && std::get<BlockName>(sp->get("k")) == "leaf");

        auto lst = std::get<std::vector<std::shared_ptr<DBNode>>>(n.get("list"));
        assert(lst.size() == 2);
        assert(std::get<BlockName>(lst[0]->get("")) == "A");
        assert(std::get<BlockName>(lst[1]->get("")) == "B");
    }

    {
        DBNode n;
        n.set(123, "lvl1", "lvl2", "value");
        assert(std::get<int>(n.get("lvl1", "lvl2", "value")) == 123);

        n.set(456, "lvl1", "lvl2", "value");
        assert(std::get<int>(n.get("lvl1", "lvl2", "value")) == 456);
    }

    {
        DBNode n;
        n.set("a", "k1");
        n.set(1, "k2");
        n.set(false, "k3");
        auto keys = n.get_keys();
        std::sort(keys.begin(), keys.end());
        assert(keys.size() == 3);
        assert(n.contains("k1"));
        assert(n.contains("k2"));
        assert(n.contains("k3"));
        assert(!n.contains("missing"));
        assert(n.countChildren() == 3);
    }

    {
        DBNode node;
        std::map<BlockName, DBNode::Value> group = {
            {BlockName("a"), 1},
            {BlockName("b"), BlockName("x")},
            {BlockName("c"), 2.2}
        };
        node.setGroup({BlockName("group1"), BlockName("sub")}, group);

        auto g = node.getGroup({BlockName("group1"), BlockName("sub")});
        assert(std::get<int>(g[BlockName("a")]) == 1);
        assert(std::get<BlockName>(g[BlockName("b")]) == "x");
        assert(std::abs(std::get<double>(g[BlockName("c")]) - 2.2) < 1e-6);

        std::map<BlockName, DBNode::Value> group2 = { {BlockName("z"), 99} };
        node.setGroup({BlockName("group1"), BlockName("sub")}, group2);
        auto g2 = node.getGroup({BlockName("group1"), BlockName("sub")});
        assert(g2.size() == 1);
        assert(std::get<int>(g2[BlockName("z")]) == 99);
    }

    {
        DBNode n;
        std::vector<std::shared_ptr<DBNode>> items;
        auto o1 = std::make_shared<DBNode>(); o1->set("X", "");
        auto o2 = std::make_shared<DBNode>(); o2->set("Y", "");
        items.push_back(o1); items.push_back(o2);
        n.set(items, "observables");

        std::stringstream ss;
        n.printJSONToStream(ss); 
        auto out = ss.str();
        assert(!out.empty());
        assert(contains_all(out, { "\"observables\"", "X", "Y" }));

    }

    {
        DBNode n;
        n.set(42, "valid");

        bool threw = false;
        try {
            (void)n.get("not_here");
        } catch (...) { threw = true; }
        assert(threw);

        threw = false;
        try {
            n.set(1, "a");
            (void)n.get("a","b");
        } catch (...) { threw = true; }
        assert(threw);
    }

    {
        DBNode n;
        n.set("streaming test", "info");
        n.set(123, "version");

        std::stringstream ss;
        n.printJSONToStream(ss);
        const auto s = ss.str();
        assert(contains_all(s, {"\"info\"", "\"streaming test\"", "\"version\"", "123"}));
    }

    std::cout << "\n All Node unit tests passed!\n";
    return 0;
}

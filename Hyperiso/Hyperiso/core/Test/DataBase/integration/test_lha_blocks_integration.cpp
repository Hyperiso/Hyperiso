#include <cassert>
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include <sstream>

#include "lha_blocks.h"
#include "LhaBlockPrototype.h"
#include "lha_elements.h"
#include "DBNode.h"

static double as_number_from_value(const DBNode::Value& v) {
    if (std::holds_alternative<double>(v)) return std::get<double>(v);
    if (std::holds_alternative<BlockName>(v)) return std::stod(std::get<BlockName>(v));
    if (std::holds_alternative<std::shared_ptr<DBNode>>(v)) {
        auto sub = std::get<std::shared_ptr<DBNode>>(v);
        auto cv = sub->get("central_value");
        if (std::holds_alternative<double>(cv)) return std::get<double>(cv);
        if (std::holds_alternative<BlockName>(cv)) return std::stod(std::get<BlockName>(cv));
    }
    throw std::runtime_error("not a numeric-compatible value");
}

int main() {
    std::cout << "== Running INTEGRATION tests for LhaBlock ==\n";

    {
        LhaBlock blk(FMASS); // itemCount=4, valueIdx=1, scaleIdx=2, rgIdx=3
        blk.readData({
            {"25", "125.0", "91.1876", "1"},
            {"35", "200.0", "91.1876", "0"}
        });

        assert(blk.get(LhaID({25})) != nullptr);
        assert(blk.get(LhaID({35})) != nullptr);

        auto n = blk.toDBNode();
        assert(n && !n->contains("scale"));
        auto g = n->getGroup({"FMASS"});
        assert(g.size() == 2);

        int seen = 0;
        for (auto& [k, v] : g) {
            double x = as_number_from_value(v);
            if (std::abs(x - 125.0) < 1e-9) ++seen;
            if (std::abs(x - 200.0) < 1e-9) ++seen;
        }
        assert(seen == 2);
    }

    {
        LhaBlock blk(GAUGE);
        blk.readData({
            {"1000.0", "1", "0.3573"},
            {"1000.0", "2", "0.6464"},
            {"1000.0", "3", "1.217"}
        });

        auto n = blk.toDBNode();
        assert(n && n->contains("scale"));
        assert(std::abs(std::get<double>(n->get("scale")) - 1000.0) < 1e-12);

        auto g = n->getGroup({"GAUGE"});
        assert(g.size() == 3);
        int ok = 0;
        for (auto& [k, v] : g) {
            double x = as_number_from_value(v);
            if (std::abs(x - 0.3573) < 1e-6) ++ok;
            if (std::abs(x - 0.6464) < 1e-6) ++ok;
            if (std::abs(x - 1.217)  < 1e-6) ++ok;
        }
        assert(ok == 3);
    }

    {
        LhaBlock blk(FCINFO);
        blk.readData({
            {"1", "GeneratorX"},
            {"2", "1.2.3"}
        });

        auto n = blk.toDBNode();
        auto g = n->getGroup({"FCINFO"});
        assert(g.size() == 2);
        int ok = 0;
        for (auto& [k, v] : g) {
            if (std::holds_alternative<std::shared_ptr<DBNode>>(v)) {
                auto sub = std::get<std::shared_ptr<DBNode>>(v);
                auto cv = sub->get("central_value");
                if (std::holds_alternative<BlockName>(cv)) {
                    auto s = std::get<BlockName>(cv);
                    if (s == "GeneratorX" || s == "1.2.3") ++ok;
                }
            }
        }
        assert(ok == 2);
    }

    std::cout << "\n All LhaBlock integration tests passed!\n";
    return 0;
}

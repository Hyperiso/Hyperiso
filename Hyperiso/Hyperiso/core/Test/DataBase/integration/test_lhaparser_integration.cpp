#include "LhaParser.h"
#include "LhaBlockPrototype.h"
#include "DBNode.h"
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <cmath>

namespace fs = std::filesystem;

static double as_num(const Node::Value& v) {
    if (std::holds_alternative<double>(v)) return std::get<double>(v);
    if (std::holds_alternative<BlockName>(v)) return std::stod(std::get<BlockName>(v));
    if (std::holds_alternative<std::shared_ptr<Node>>(v)) {
        auto sub = std::get<std::shared_ptr<Node>>(v);
        auto cv = sub->get("central_value");
        if (std::holds_alternative<double>(cv)) return std::get<double>(cv);
        if (std::holds_alternative<BlockName>(cv)) return std::stod(std::get<BlockName>(cv));
    }
    throw std::runtime_error("not numeric");
}

static void create_flha(const fs::path& path) {
    std::ofstream f(path);
    f <<
        "Block GAUGE  Q= 1.00000000E+03\n"
        "  1  3.57300000E-01  # g1\n"
        "  2  6.46400000E-01  # g2\n"
        "  3  1.21700000E+00  # g3\n"
        "BLOCK MASS\n"
        "  25  1.25000000E+02\n"
        "BLOCK FMASS\n"
        "  13  1.234  91.1876  1\n";
}

int main() {
    std::cout << "== Running INTEGRATION tests for LhaParser ==\n";

    const fs::path file = "parser_it.flha";
    create_flha(file);

    LhaParser p;
    p.set_prototypes({GAUGE, MASS, FMASS});

    auto root = p.readFromFile(file.string());

    // GAUGE
    auto GA = std::get<std::shared_ptr<Node>>(root->get("GAUGE"));
    assert(std::abs(as_num(GA->get("1")) - 0.3573) < 1e-4);
    assert(std::abs(as_num(GA->get("2")) - 0.6464) < 1e-4);
    assert(std::abs(as_num(GA->get("3")) - 1.217)  < 1e-6);
    assert(std::abs(std::get<double>(GA->get("scale")) - 1000.0) < 1e-9);

    // MASS
    auto MA = std::get<std::shared_ptr<Node>>(root->get("MASS"));
    assert(std::abs(as_num(MA->get("25")) - 125.0) < 1e-2);

    // FMASS (scale par-élément)
    auto FM = std::get<std::shared_ptr<Node>>(root->get("FMASS"));
    auto mu = std::get<std::shared_ptr<Node>>(FM->get("13"));
    assert(std::abs(std::get<double>(mu->get("central_value")) - 1.234) < 1e-12);
    assert(std::abs(std::get<double>(mu->get("scale")) - 91.1876) < 1e-6);
    assert(std::get<int>(mu->get("renormalization_scheme")) == 1);

    fs::remove(file);

    std::cout << "\n✅ LhaParser integration tests passed!\n";
    return 0;
}

#include "LhaParser.h"
#include "LhaBlockPrototype.h"
#include "DBNode.h"
#include <cassert>
#include <iostream>
#include <cmath>
#include <sstream>
#include <unordered_set>

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

int main() {
    std::cout << "== Running UNIT tests for LhaParser ==\n";

    // Prépare un parser avec les prototypes utiles
    LhaParser p;
    std::unordered_set<Prototype> protos = {
        GAUGE, MASS, FMASS  // GAUGE globalScale=true ; MASS simple ; FMASS avec scaleIdx/rgIdx
    };
    p.set_prototypes(protos);

    // Source LHA minimaliste
    const std::string src =
        "Block  GAUGE   Q= 1000\n"
        "  1   3.57300000E-01   # g1\n"
        "  3   1.21700000E+00   # g3\n"
        "BLOCK MASS\n"
        "  25  1.25000000E+02\n"
        "DeCaY 25  1.0\n"            // doit être ignoré
        "Block FMASS\n"
        "  13  1.234  91.1876  1\n"  // value=1.234, scale=91.1876, rg=1
        "# end\n";

    auto root = p.parse(src);

    // 1) GAUGE : présent + scale global stocké sous GAUGE/scale
    assert(root->contains("GAUGE"));
    auto G = std::get<std::shared_ptr<Node>>(root->get("GAUGE"));
    auto g1 = G->get("1");
    auto g3 = G->get("3");
    assert(std::abs(as_num(g1) - 0.3573) < 1e-4);
    assert(std::abs(as_num(g3) - 1.217)  < 1e-6);
    assert(std::holds_alternative<double>(G->get("scale")));
    assert(std::abs(std::get<double>(G->get("scale")) - 1000.0) < 1e-12);

    // 2) MASS : présent, pas de scale global, valeur OK
    assert(root->contains("MASS"));
    auto M = std::get<std::shared_ptr<Node>>(root->get("MASS"));
    auto h = M->get("25");
    assert(std::abs(as_num(h) - 125.0) < 1e-2);

    // 3) FMASS : valeur + scale par élément + rg par élément
    assert(root->contains("FMASS"));
    auto F = std::get<std::shared_ptr<Node>>(root->get("FMASS"));
    auto mu = std::get<std::shared_ptr<Node>>(F->get("13")); // sous-noeud
    assert(std::abs(std::get<double>(mu->get("central_value")) - 1.234) < 1e-12);
    assert(std::abs(std::get<double>(mu->get("scale")) - 91.1876) < 1e-6);
    assert(std::get<int>(mu->get("renormalization_scheme")) == 1);

    std::cout << "\n✅ LhaParser unit tests passed!\n";
    return 0;
}

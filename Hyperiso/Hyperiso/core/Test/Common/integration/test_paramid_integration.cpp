#include <cassert>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <algorithm>

// Adapte ces includes à ton arborescence :
#include "General.h"


int main() {
    std::cout << "== Running INTEGRATION tests for ParamId ==\n";

    auto Tmass     = static_cast<ParameterType>(10);
    auto Tcoupling = static_cast<ParameterType>(20);

    // 1) Utilisation comme clé d’une table triée (std::map)
    {
        std::map<ParamId, double> table;

        table.emplace(ParamId(Tmass,     BlockName("MASS"),  LhaID(25)), 125.0);
        table.emplace(ParamId(Tmass,     BlockName("MASS"),  LhaID(35)), 200.0);
        table.emplace(ParamId(Tcoupling, BlockName("GAUGE"), LhaID(1)),  0.3573);

        // Accès
        assert(table.find(ParamId(Tmass, BlockName("MASS"), LhaID(25))) != table.end());
        assert(table.at(ParamId(Tmass, BlockName("MASS"), LhaID(35))) == 200.0);

        // Tri (par type, block, code) — on vérifie les extrémités
        auto first = table.begin()->first;
        auto last  = std::prev(table.end())->first;

        // nullopt < Tmass < Tcoupling ; ici pas de nullopt, donc Tmass < Tcoupling
        assert(first.type && *first.type == Tmass);
        assert(last.type  && *last.type  == Tcoupling);
    }

    // 2) Utilisation en unordered_map (OK tant que BlockName n’a pas d’alias “bizarres”)
    {
        std::unordered_map<ParamId, const char*> cache;
        ParamId p1(Tmass, BlockName("MASS"),  LhaID(25));
        ParamId p2(Tcoupling, BlockName("GAUGE"), LhaID(3));

        cache[p1] = "mh";
        cache[p2] = "g1";

        assert(std::string(cache[p1]) == "mh");
        assert(std::string(cache[p2]) == "g1");

        // Remplacement
        cache[p1] = "mh_updated";
        assert(std::string(cache[p1]) == "mh_updated");
    }

    std::cout << "\n✅ ParamId integration tests passed!\n";
    return 0;
}

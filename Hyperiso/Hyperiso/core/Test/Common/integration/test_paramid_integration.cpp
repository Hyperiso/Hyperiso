#include <cassert>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <algorithm>

#include "Include.h"


int main() {
    std::cout << "== Running INTEGRATION tests for ParamId ==\n";

    auto Tmass     = static_cast<ParameterType>(10);
    auto Tcoupling = static_cast<ParameterType>(20);

    {
        std::map<ParamId, double> table;

        table.emplace(ParamId(Tmass,     BlockName("MASS"),  LhaID(25)), 125.0);
        table.emplace(ParamId(Tmass,     BlockName("MASS"),  LhaID(35)), 200.0);
        table.emplace(ParamId(Tcoupling, BlockName("GAUGE"), LhaID(1)),  0.3573);

        assert(table.find(ParamId(Tmass, BlockName("MASS"), LhaID(25))) != table.end());
        assert(table.at(ParamId(Tmass, BlockName("MASS"), LhaID(35))) == 200.0);

        auto first = table.begin()->first;
        auto last  = std::prev(table.end())->first;

        assert(first.type && *first.type == Tmass);
        assert(last.type  && *last.type  == Tcoupling);
    }

    {
        std::unordered_map<ParamId, const char*> cache;
        ParamId p1(Tmass, BlockName("MASS"),  LhaID(25));
        ParamId p2(Tcoupling, BlockName("GAUGE"), LhaID(3));

        cache[p1] = "mh";
        cache[p2] = "g1";

        assert(std::string(cache[p1]) == "mh");
        assert(std::string(cache[p2]) == "g1");

        cache[p1] = "mh_updated";
        assert(std::string(cache[p1]) == "mh_updated");
    }

    std::cout << "\n ParamId integration tests passed!\n";
    return 0;
}

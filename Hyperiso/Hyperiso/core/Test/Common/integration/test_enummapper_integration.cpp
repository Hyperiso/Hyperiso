#include <cassert>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

#include "EnumMapper.h"

// Reprend la même définition de Fruit / FruitMap que dans l’unit (ou mets-les dans un header commun).
enum class Fruit { Apple = 0, Banana = 1, Cherry = 2 };

struct FruitMap : EnumMapperBase<Fruit, FruitMap> {
    static const std::map<Fruit, std::string>& mapping() {
        static const std::map<Fruit, std::string> m = {
            {Fruit::Apple,  "APPLE"},
            {Fruit::Banana, "BANANA"},
            {Fruit::Cherry, "CHERRY"},
        };
        return m;
    }
    static const std::map<std::string, Fruit>& inverse_mapping() {
        static const std::map<std::string, Fruit> inv = invert_map(mapping());
        return inv;
    }
};

static bool contains_all(const std::vector<std::string>& v, std::initializer_list<const char*> need) {
    for (auto* s : need) if (std::find(v.begin(), v.end(), std::string(s)) == v.end()) return false;
    return true;
}

int main() {
    std::cout << "== Running INTEGRATION tests for EnumMapperBase ==\n";

    // 1) Round-trip vector<Enum> -> vector<string> -> vector<Enum>
    {
        std::vector<Fruit> in = {Fruit::Banana, Fruit::Apple, Fruit::Cherry};

        std::vector<std::string> names;
        names.reserve(in.size());
        for (auto e : in) names.push_back(FruitMap::str(e));

        assert(contains_all(names, {"BANANA","APPLE","CHERRY"}));

        std::vector<Fruit> out;
        out.reserve(names.size());
        for (auto& s : names) out.push_back(FruitMap::enum_elt(s));

        assert(in.size() == out.size());
        for (size_t i=0;i<in.size();++i) assert(in[i] == out[i]);
    }

    // 2) Intégration dans une map typée par l’enum
    {
        std::map<Fruit, double> prices = {
            {Fruit::Apple,  1.0},
            {Fruit::Banana, 0.8},
            {Fruit::Cherry, 2.5},
        };

        // affichage textuel via mapper + vérifs
        auto nameA = FruitMap::str(Fruit::Apple);
        auto nameB = FruitMap::str(Fruit::Banana);
        auto nameC = FruitMap::str(Fruit::Cherry);
        assert(nameA == "APPLE" && nameB == "BANANA" && nameC == "CHERRY");

        // lookup côté inverse
        assert(FruitMap::enum_elt(nameA) == Fruit::Apple);
        assert(FruitMap::enum_elt(nameB) == Fruit::Banana);
        assert(FruitMap::enum_elt(nameC) == Fruit::Cherry);

        // vérif rapide
        assert(prices.at(Fruit::Cherry) > prices.at(Fruit::Apple));
    }

    std::cout << "\n✅ EnumMapperBase integration tests passed!\n";
    return 0;
}

#include <cassert>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <stdexcept>

#include "EnumMapper.h"

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

template <typename T>
static bool vec_eq(const std::vector<T>& a, const std::vector<T>& b) {
    return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin());
}

int main() {
    std::cout << "== Running UNIT tests for EnumMapperBase ==\n";

    {
        assert(FruitMap::str(Fruit::Apple)  == "APPLE");
        assert(FruitMap::str(Fruit::Banana) == "BANANA");
        assert(FruitMap::str(Fruit::Cherry) == "CHERRY");
    }

    {
        assert(FruitMap::enum_elt("APPLE")  == Fruit::Apple);
        assert(FruitMap::enum_elt("BANANA") == Fruit::Banana);
        assert(FruitMap::enum_elt("CHERRY") == Fruit::Cherry);

        bool threw = false;
        try {
            (void)FruitMap::enum_elt("PEAR");
        } catch (const std::out_of_range&) { threw = true; }
        assert(threw);
    }

    {
        auto s = FruitMap::get_str();
        std::vector<std::string> expected = {"APPLE","BANANA","CHERRY"};
        assert(vec_eq(s, expected));
    }

    {
        auto e = FruitMap::get_enum();
        std::vector<Fruit> expected = {Fruit::Apple, Fruit::Banana, Fruit::Cherry};
        assert(vec_eq(e, expected));
    }

    {
        std::map<int, std::string> orig = {{1,"one"},{2,"two"},{3,"three"}};
        auto inv = invert_map(orig);
        assert(inv.at("one")   == 1);
        assert(inv.at("two")   == 2);
        assert(inv.at("three") == 3);
    }

    {
        std::map<int, std::string> orig = {{1,"x"},{2,"x"}};
        auto inv = invert_map(orig);
        assert(inv.size() == 1);
        assert(inv.at("x") == 2);
    }

    std::cout << "\n All EnumMapperBase unit tests passed!\n";
    return 0;
}

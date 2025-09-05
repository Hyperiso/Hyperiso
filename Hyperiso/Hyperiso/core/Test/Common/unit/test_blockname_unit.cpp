// test_blockname_unit.cpp
#include <cassert>
#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <sstream>

// adapte l'include si besoin
#include "General.h"

static bool in_set(const std::unordered_set<std::string>& S, const std::string& x) {
    return S.find(x) != S.end();
}

static std::unordered_set<std::string> split_slash(const std::string& s) {
    std::unordered_set<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == '/') { if (!cur.empty()) { out.insert(cur); cur.clear(); } }
        else cur.push_back(c);
    }
    if (!cur.empty()) out.insert(cur);
    return out;
}

int main() {
    std::cout << "== Running UNIT tests for BlockName ==\n";

    // 1) Constructeurs
    {
        BlockName empty;
        assert(std::string(empty) == ""); // to_string() via operator std::string()

        BlockName s1(std::string("MASS"));
        assert(s1 == "MASS");

        BlockName s2("GAUGE");
        assert(s2 == std::string("GAUGE"));

        BlockName many{ "YU", "UCOUPL" };
        assert(many == "YU");
        assert(many == "UCOUPL");
        assert(many != "YD");

        std::unordered_set<std::string> aliases = {"A","B","C"};
        BlockName from_set(aliases);
        for (auto& a : aliases) assert(from_set == a);
    }

    // 2) to_upper()
    {
        BlockName b{"mass","Gauge"};
        b.to_upper();
        assert(b == "MASS");
        assert(b == "GAUGE");
        assert(b != "mass");
    }

    // 3) hasAlias / addAlias
    {
        BlockName b("MASS");
        assert(b.hasAlias("MASS"));
        assert(!b.hasAlias("GAUGE"));
        b.addAlias("HMM");
        assert(b.hasAlias("HMM"));
    }

    // 4) operator== / operator!= (avec BlockName et std::string)
    {
        BlockName a{"YU","UCOUPL"};
        BlockName b{"UCOUPL","XYZ"};
        BlockName c{"YD"};

        // égal si intersection non-vide
        assert(a == b);
        assert(!(a == c));
        assert(a != c);

        assert(a == std::string("YU"));
        assert(std::string("UCOUPL") == a);
        assert(a != std::string("YD"));
    }

    // 5) operator+ (préfixe string)
    {
        BlockName b{"MASS","GAUGE"};
        BlockName pref = std::string("SUSY_") + b;
        auto al = pref.get_alias();
        assert(in_set(al, "SUSY_MASS"));
        assert(in_set(al, "SUSY_GAUGE"));
    }

    // 6) operator<< (ordre des alias non garanti → on vérifie l’ensemble)
    {
        BlockName b{"NMIX", "UMIX", "VMIX"};
        std::stringstream ss;
        ss << b;
        auto tokens = split_slash(ss.str());
        assert(tokens.size() == 3);
        assert(in_set(tokens, "NMIX"));
        assert(in_set(tokens, "UMIX"));
        assert(in_set(tokens, "VMIX"));
    }

    // 7) operator< (tri déterministe via sets)
    {
        std::vector<BlockName> v = {
            BlockName{"B"}, BlockName{"A"}, BlockName{"A","A2"}, BlockName{"C"}
        };
        std::sort(v.begin(), v.end());
        // ordre attendu (lexico sur sets) : {"A"}, {"A","A2"}, {"B"}, {"C"}
        assert(v[0] == "A");
        assert(v[1] == "A2" || v[1] == "A"); // contient A2 ET A; l'égalité renvoie true si intersection
        assert(v[2] == "B");
        assert(v[3] == "C");
    }

    // 8) conversion implicite -> std::string : doit renvoyer une des aliases
    {
        BlockName b{"MASS","GAUGE"};
        std::string one = b; // operator std::string()
        auto al = b.get_alias();
        assert(in_set(al, one)); // c'est l'une des entrées (ordre non spécifié)
    }

    std::cout << "\n✅ All BlockName unit tests passed!\n";
    return 0;
}

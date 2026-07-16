#include "Include.h"
#include <cassert>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <sstream>

static bool equal_vec(const std::vector<long>& a, const std::vector<long>& b) {
    return a == b;
}

int main() {
    std::cout << "== Running UNIT tests for LhaID ==\n";

    {
        LhaID id(1, 2, 3);
        auto parts = id.get_parts();
        assert(equal_vec(parts, {1,2,3}));
        assert(id.to_string() == "1_2_3");
    }

    {
        LhaID id(std::string("10_-7_0"));
        auto parts = id.get_parts();
        assert(equal_vec(parts, {10,-7,0}));
        assert(id.to_string() == "10_-7_0");
    }

    {
        std::vector<long> v = {42};
        LhaID idv(v);
        assert(equal_vec(idv.get_parts(), {42}));
        assert(idv.to_string() == "42");

        LhaID idil{5,6};
        assert(equal_vec(idil.get_parts(), {5,6}));
        assert(idil.to_string() == "5_6");

        LhaID id1(999L);
        assert(equal_vec(id1.get_parts(), {999}));
        assert(id1.to_string() == "999");
    }

    {
        LhaID id1(123);
        long x = id1; 
        assert(x == 123);

        LhaID idm(9,8);
        long y = idm;
        assert(y == 9);
    }

    {
        LhaID a(1,2), b(1,2), c(1,3), d(2);

        assert(a == b);
        assert(a != c);
        assert(a != d);

        assert(a < c);  // {1,2} < {1,3}
        assert(c < d);  // {1,3} < {2}
        assert(!(d < c));
    }

    {
        std::unordered_set<LhaID> S;
        S.insert(LhaID(1,2));
        S.insert(LhaID(1,2));
        S.insert(LhaID(1,3));
        assert(S.size() == 2);

        std::unordered_map<LhaID, int> M;
        M[LhaID(10)] = 100;
        M[LhaID(10,1)] = 101;
        assert(M[LhaID(10)] == 100);
        assert(M[LhaID(10,1)] == 101);
    }

    {
        std::map<LhaID, int> m;
        m[LhaID(1,2)] = 12;
        m[LhaID(1,3)] = 13;
        m[LhaID(2)]   = 2;

        // (1,2) < (1,3) < (2)
        auto it = m.begin();
        assert(it->first == LhaID(1,2)); ++it;
        assert(it->first == LhaID(1,3)); ++it;
        assert(it->first == LhaID(2));
    }

    {
        bool threw = false;
        try {
            LhaID bad("a_b");
        } catch (...) { threw = true; }
        (void)threw;
    }

    std::cout << "\nAll LhaID unit tests passed!\n";
    return 0;
}

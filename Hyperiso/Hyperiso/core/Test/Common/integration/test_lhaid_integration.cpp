#include "Include.h"
#include <cassert>
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <map>

static std::string join(const std::vector<LhaID>& v) {
    std::string s;
    for (size_t i=0;i<v.size();++i) {
        s += v[i].to_string();
        if (i+1<v.size()) s += ",";
    }
    return s;
}

int main() {
    std::cout << "== Running INTEGRATION tests for LhaID ==\n";

    {
        std::vector<LhaID> ids = {
            LhaID(2),
            LhaID(1,3),
            LhaID(std::string("1_2")),
            LhaID{1,1,1},
            LhaID(1,2,3)
        };
        std::sort(ids.begin(), ids.end());
        auto seq = join(ids);
        auto pos = seq.find("1_1_1");
        assert(pos != std::string::npos);
        pos = seq.find("1_2", pos);
        assert(pos != std::string::npos);
        pos = seq.find("1_2_3", pos);
        assert(pos != std::string::npos);
        pos = seq.find("1_3", pos);
        assert(pos != std::string::npos);
        pos = seq.find("2", pos);
        assert(pos != std::string::npos);
    }

    {
        std::unordered_map<LhaID, const char*> U;
        U[LhaID(1,2)] = "A";
        U[LhaID(1,3)] = "B";
        U[LhaID(1,2)] = "A2"; 
        assert(std::string(U[LhaID(1,2)]) == "A2");
        assert(std::string(U[LhaID(1,3)]) == "B");

        std::map<LhaID, int> M;
        M[LhaID(1,2)] = 12;
        M[LhaID(1,3)] = 13;
        M[LhaID(2)]   = 2;
        auto it = M.begin();
        assert(it->first == LhaID(1,2)); ++it;
        assert(it->first == LhaID(1,3)); ++it;
        assert(it->first == LhaID(2));
    }

    std::cout << "\n LhaID integration tests passed!\n";
    return 0;
}

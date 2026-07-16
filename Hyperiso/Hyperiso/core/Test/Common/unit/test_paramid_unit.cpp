#include <cassert>
#include <iostream>
#include <unordered_map>
#include <map>
#include <vector>
#include <algorithm>

#include "Include.h"

static bool is_nullopt(const std::optional<ParameterType>& t) {
    return !t.has_value();
}

int main() {
    std::cout << "== Running UNIT tests for ParamId ==\n";

    {
        ParamId pid;
        assert(is_nullopt(pid.type));
        assert(static_cast<std::string>(pid.block) == "NULL");
        assert(static_cast<long>(pid.code) == 0);
    }

    {
        ParamId pid(BlockName("MASS"), LhaID(25));
        assert(is_nullopt(pid.type));
        assert(pid.block == "MASS");
        assert(pid.code == LhaID(25));
    }

    {
        auto T1 = static_cast<ParameterType>(1);
        auto T2 = static_cast<ParameterType>(2);

        ParamId pid(T1, BlockName("GAUGE"), LhaID(3));
        assert(pid.type.has_value() && *pid.type == T1);
        assert(pid.block == "GAUGE");
        assert(pid.code  == LhaID(3));

        pid.set_parameter_type(T2);
        assert(pid.type.has_value() && *pid.type == T2);
    }

    {
        auto T1 = static_cast<ParameterType>(1);

        ParamId a(T1, BlockName("MASS"),  LhaID(25));
        ParamId b(T1, BlockName("MASS"),  LhaID(25));
        ParamId c(T1, BlockName("MASS"),  LhaID(35));
        ParamId d(    BlockName("MASS"),  LhaID(25));
        ParamId e(T1, BlockName("GAUGE"), LhaID(25));

        assert(a == b);
        assert(a != c);
        assert(a != d);
        assert(a != e);
    }

    {
        auto T1 = static_cast<ParameterType>(1);
        auto T2 = static_cast<ParameterType>(2);

        std::vector<ParamId> v = {
            ParamId(        BlockName("GAUGE"), LhaID(1)),
            ParamId(T2,     BlockName("GAUGE"), LhaID(1)),
            ParamId(T1,     BlockName("GAUGE"), LhaID(2)),
            ParamId(T1,     BlockName("GAUGE"), LhaID(1)),
            ParamId(T1,     BlockName("MASS"),  LhaID(1))
        };
        std::sort(v.begin(), v.end());

        //  - (nullopt, GAUGE, 1)
        //  - (T1, GAUGE, 1)
        //  - (T1, GAUGE, 2)
        //  - (T1, MASS,  1)
        //  - (T2, GAUGE, 1)
        assert(is_nullopt(v[0].type) && v[0].block == "GAUGE" && v[0].code == LhaID(1));
        assert(v[1].type && *v[1].type == T1 && v[1].block == "GAUGE" && v[1].code == LhaID(1));
        assert(v[2].type && *v[2].type == T1 && v[2].block == "GAUGE" && v[2].code == LhaID(2));
        //  "GAUGE" < "MASS"
        assert(v[3].type && *v[3].type == T1 && v[3].block == "MASS"  && v[3].code == LhaID(1));
        assert(v[4].type && *v[4].type == T2 && v[4].block == "GAUGE" && v[4].code == LhaID(1));
    }

    {
        auto T1 = static_cast<ParameterType>(1);
        ParamId k1(T1, BlockName("MASS"),  LhaID(25));
        ParamId k2(T1, BlockName("GAUGE"), LhaID(3));

        std::unordered_map<ParamId, int> M;
        M[k1] = 125;
        M[k2] = 42;

        assert(M[k1] == 125);
        assert(M[k2] == 42);
    }

    std::cout << "\nAll ParamId unit tests passed!\n";
    return 0;
}

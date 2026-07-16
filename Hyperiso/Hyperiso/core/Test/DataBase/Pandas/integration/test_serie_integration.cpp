#include <cassert>
#include <iostream>
#include <vector>
#include "Series.h"

int main() {
    std::cout << "== Series INTEGRATION ==\n";

    Series<int>    si("ints");
    Series<double> sd("dbls");
    Series<std::string> ss("strs");

    for (int i = 1; i <= 5; ++i) {
        si.add(i);
        sd.add(0.5 * i);
        ss.add(std::to_string(i));
    }

    assert(si.size() == 5);
    assert(sd.size() == 5);
    assert(ss.size() == 5);

    assert(si.min() == 1 && si.max() == 5);
    assert(sd.min() == 0.5 && sd.max() == 2.5);

    Series<double> sx("with_index");
    sx.add(1.1, "k1");
    sx.add(2.2, "k2");
    sx.add(3.3, "k3");
    assert(std::abs(sx.at("k2") - 2.2) < 1e-12);

    std::cout << "✅ INTEGRATION OK\n";
    return 0;
}

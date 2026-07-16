#include <cassert>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include "Series.h"

int main() {
    std::cout << "== Series UNIT ==\n";

    Series<double> s("S");
    s.add(1.0);
    s.add(2.0);
    s.add(3.0);
    assert(s.size() == 3);
    assert(std::abs(s.iat(0) - 1.0) < 1e-12);
    assert(std::abs(s.iat(2) - 3.0) < 1e-12);

    bool threw = false;
    try { (void)s.iat(3); } catch (const std::out_of_range&) { threw = true; }
    assert(threw);

    assert(std::abs(s.min() - 1.0) < 1e-12);
    assert(std::abs(s.max() - 3.0) < 1e-12);
    assert(std::abs(s.mean() - 2.0) < 1e-12);
    assert(std::abs(s.stddev() - std::sqrt((1+0+1)/3.0)) < 1e-12);

    auto qs = s.quartiles();
    assert(qs.size() == 3);
    assert(std::abs(qs[0] - 1.0) < 1e-12);
    assert(std::abs(qs[1] - 2.0) < 1e-12);
    assert(std::abs(qs[2] - 3.0) < 1e-12);

    Series<int> si("I");
    si.add(10, "a");
    si.add(20, "b");
    si.add(30, "c");
    assert(si.size() == 3);
    assert(si.at("a") == 10);
    assert(si.at("c") == 30);
    bool threw_idx = false;
    try { (void)si.at("z"); } catch (const std::invalid_argument&) { threw_idx = true; }
    assert(threw_idx);

    std::cout << "✅ UNIT OK\n";
    return 0;
}

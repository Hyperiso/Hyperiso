#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <vector>

#include "MemoryManager.h"
#include "WilsonInterface.h"
#include "BlockProxy.h"
#include "wcoef_ids.hpp"

namespace fs = std::filesystem;

static bool approx(double a, double b, double eps = 1e-6) {
    return std::fabs(a - b) <= eps;
}

static fs::path write_temp_fwcoef_tot_only() {
    fs::path p = fs::temp_directory_path() / "fwcoef_tot_only_test.flha";
    std::ofstream out(p);
    out <<
R"(Block FWCOEF Q= 1.60846e+02 # Wilson coefficients at scale Q
#id                 order   M   value           comment
03040405    6161    00      2   0.00000000e+00  # C1^0
03040405    4141    00      2   1.00000000e+00  # C2^0
0305        4422    00      2   -1.53496321e-01 # C7^0
0305        6421    00      2   -9.51462419e-02 # C8^0
)";
    return p;
}

int main() {
    std::cout << "== FWCOEF integration (TOT-only) ==\n";

    auto hyp = HyperisoMaster();
    HyperisoConfig cfg;
    cfg.flags[ExternalFlag::HAS_WILSON_INPUT] = true;

    fs::path flha = write_temp_fwcoef_tot_only();
    hyp.init(flha.string(), cfg);

    WilsonInterface wi;
    WilsonBuildConfig build_cfg({WGroup::B}, 81.0, 4.18, QCDOrder::LO);
    wi.build(build_cfg);

    auto idC7 = WCoefMapper::enum_of(WCoefMapper::enum_elt("C7")).value();
    auto idC8 = WCoefMapper::enum_of(WCoefMapper::enum_elt("C8")).value();

    complex_t C7_tot = wi.getMatchingCoefficient(WGroup::B, idC7, QCDOrder::LO, ContributionType::TOTAL);
    complex_t C8_tot = wi.getMatchingCoefficient(WGroup::B, idC8, QCDOrder::LO, ContributionType::TOTAL);

    assert(approx(C7_tot.real(), -1.53496321e-01));
    assert(approx(C7_tot.imag(), 0.0));
    assert(approx(C8_tot.real(), -9.51462419e-02));
    assert(approx(C8_tot.imag(), 0.0));

    complex_t C7_sm  = wi.getMatchingCoefficient(WGroup::B, idC7, QCDOrder::LO, ContributionType::SM);
    complex_t C7_bsm = wi.getMatchingCoefficient(WGroup::B, idC7, QCDOrder::LO, ContributionType::BSM);

    assert(approx((C7_sm + C7_bsm).real(), C7_tot.real(), 1e-6));
    assert(approx((C7_sm + C7_bsm).imag(), C7_tot.imag(), 1e-6));

    auto idC9 = WCoefMapper::enum_of(WCoefMapper::enum_elt("C9")).value();

    complex_t C9_sm  = wi.getMatchingCoefficient(WGroup::B, idC9, QCDOrder::LO, ContributionType::SM);
    complex_t C9_bsm = wi.getMatchingCoefficient(WGroup::B, idC9, QCDOrder::LO, ContributionType::BSM);
    complex_t C9_tot = wi.getMatchingCoefficient(WGroup::B, idC9, QCDOrder::LO, ContributionType::TOTAL);

    assert(approx(C9_bsm.real(), 0.0));
    assert(approx(C9_bsm.imag(), 0.0));
    assert(approx(C9_tot.real(), C9_sm.real(), 1e-6));
    assert(approx(C9_tot.imag(), C9_sm.imag(), 1e-6));

    std::cout << "FWCOEF TOT-only integration suite passed.\n";
    return 0;
}

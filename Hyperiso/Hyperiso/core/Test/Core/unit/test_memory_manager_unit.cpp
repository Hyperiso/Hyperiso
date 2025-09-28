#include <cassert>
#include <iostream>
#include <algorithm>

#include "HyperisoMaster.h"
#include "MemoryManager.h"
#include "Config.h"

#include "dummies_mm.hpp"

static bool hasType(const std::vector<ParameterType>& v, ParameterType t){
    return std::find(v.begin(), v.end(), t) != v.end();
}

int main() {
    std::cout << "== HyperisoMaster UNIT ==\n";

    const fs::path lha_rel = "unit_master.lha";
    fs::path sandbox = prepare_assets_for_mm(lha_rel);

    auto ba   = std::make_shared<DummyBA>();
    auto cp   = std::make_shared<DummyCorr<ParamId>>();
    auto co   = std::make_shared<DummyCorr<ObservableId>>();
    auto spec = std::make_shared<DummySpectrum>();
    auto paths= std::make_shared<TestPathsProvider>(sandbox);
    (void)MemoryManager::Create(ba, cp, co, spec, paths);

    HyperisoMaster master;
    Config cfg;
    cfg.model = Model::SM;
    cfg.flags[ExternalFlag::HAS_WILSON_INPUT] = true;
    cfg.flags[ExternalFlag::USE_MARTY] = false;

    master.init(lha_rel.string(), cfg);

    assert(master.check_flag(ExternalFlag::HAS_WILSON_INPUT) == true);
    assert(master.check_flag(ExternalFlag::USE_MARTY) == false);
    assert(master.get_model() == Model::SM);

    const auto& types = MemoryManager::GetInstance()->getMemoryCache().parameter_types;
    assert(hasType(types, ParameterType::WILSON));

    master.init(lha_rel.string());
    assert(master.get_model() == Model::SM);

    Config cfg2 = cfg;
    cfg2.flags[ExternalFlag::USE_MARTY] = true;
    master.switch_lha(lha_rel.string(), cfg2);
    assert(master.check_flag(ExternalFlag::USE_MARTY) == true);

    std::cout << " UNIT OK\n";
    return 0;
}

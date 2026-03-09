#include <cassert>
#include <iostream>
#include <algorithm>

#include "HyperisoMaster.h"
#include "MemoryManager.h"
#include "Config.h"
#include "BlockAccessor.h"

#include "dummies_mm.hpp"

static bool hasType(const std::vector<ParameterType>& v, ParameterType t){
    return std::find(v.begin(), v.end(), t) != v.end();
}

static double firstVal(const std::shared_ptr<BlockAccessor>& ba, const BlockName& b){
    auto& blk = ba->at(b);
    if (blk->getItems().empty()) return 0.0;
    return blk->getItems().begin()->second->get_val();
}

int main() {
    std::cout << "== HyperisoMaster INTEGRATION ==\n";

    const fs::path lha_rel1 = "int_master_1.lha";
    const fs::path lha_rel2 = "int_master_2.lha";
    fs::path sandbox = prepare_assets_for_mm(lha_rel1);
    prepare_assets_for_mm(lha_rel2);

    auto ba   = std::make_shared<DummyBA>();
    auto cp   = std::make_shared<DummyCorr<ParamId>>();
    auto co   = std::make_shared<DummyCorr<BinnedObservableId>>();
    auto spec = std::make_shared<DummySpectrum>();
    auto paths= std::make_shared<TestPathsProvider>(sandbox);
    (void)MemoryManager::Create(ba, cp, co, spec, paths);

    HyperisoMaster master;
    HyperisoConfig cfg;
    cfg.model = Model::SM;
    cfg.flags[ExternalFlag::HAS_WILSON_INPUT] = false;
    // cfg.flags[ExternalFlag::USE_MARTY] = false;

    master.init(lha_rel1.string(), cfg);
    assert(master.get_model() == Model::SM);

    auto* mm = MemoryManager::GetInstance();

    auto slice1 = mm->extract_blocks({"MASS","GAUGE"});
    assert(slice1->contains("MASS"));
    assert(slice1->contains("GAUGE"));

    {
        const auto& types = mm->getMemoryCache().parameter_types;
        assert(hasType(types, ParameterType::WILSON));
        assert(!hasType(types, ParameterType::BSM));
    }

    double v1 = firstVal(slice1, "GAUGE");

    HyperisoConfig cfg_sw = cfg;
    // cfg_sw.flags[ExternalFlag::USE_MARTY] = true;
    master.switch_lha(lha_rel2.string(), cfg_sw);
    // assert(master.check_flag(ExternalFlag::USE_MARTY) == true);
    //TODO : redo tests
    auto slice2 = mm->extract_blocks({"MASS","GAUGE"});
    assert(slice2->contains("MASS"));
    assert(slice2->contains("GAUGE"));
    double v2 = firstVal(slice2, "GAUGE");
    std::cout << "v1 "<< v1 << std::endl;
    std::cout << "v2 " << v2 << std::endl; 
    assert(v1 != v2);

    {
        const auto& types = mm->getMemoryCache().parameter_types;
        assert(hasType(types, ParameterType::WILSON));
        assert(!hasType(types, ParameterType::BSM));
    }

    mm->switch_model(Model::THDM);
    {
        const auto& types = mm->getMemoryCache().parameter_types;
        assert(hasType(types, ParameterType::WILSON));
        assert(hasType(types, ParameterType::BSM));
        assert(std::count(types.begin(), types.end(), ParameterType::WILSON) == 1);
        assert(std::count(types.begin(), types.end(), ParameterType::BSM) == 1);
    }

    std::cout << " INTEGRATION OK\n";
    return 0;
}

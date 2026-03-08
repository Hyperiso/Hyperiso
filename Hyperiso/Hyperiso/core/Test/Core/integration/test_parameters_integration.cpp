#include <cassert>
#include <iostream>
#include "Parameters.h"
#include "MemoryManager.h"
#include "Config.h"
#include "Paths.h"
#include "BlockAccessor.h"
#include "dummies_mm.hpp"

static bool hasType(const std::vector<ParameterType>& v, ParameterType t){
  return std::find(v.begin(), v.end(), t) != v.end();
}

int main() {
  std::cout << "== Parameters INTEGRATION ==\n";

  // const fs::path lha1 = "int_params_1.lha";
  // const fs::path lha2 = "int_params_2.lha";
  // prepare_assets_for_mm(lha1);
  // prepare_assets_for_mm(lha2);

    const fs::path lha1 = "int_params_1.lha";
    const fs::path lha2 = "int_params_2.lha";
    fs::path sandbox = prepare_assets_for_mm(lha1);
    prepare_assets_for_mm(lha2);

  auto ba   = std::make_shared<DummyBA>();
  auto cp   = std::make_shared<DummyCorr<ParamId>>();
  auto co   = std::make_shared<DummyCorr<BinnedObservableId>>();
  auto spec = std::make_shared<DummySpectrum>();
  auto tpp = std::make_shared<TestPathsProvider>(sandbox);
  (void)MemoryManager::Create(ba, cp, co, spec, tpp);

  HyperisoConfig cfg;
  cfg.model = Model::SM;
  cfg.flags[ExternalFlag::HAS_WILSON_INPUT] = false;
  // cfg.flags[ExternalFlag::USE_MARTY] = false;

  auto* mm = MemoryManager::GetInstance();
  mm->init(lha1.string(), cfg);

  auto sm1 = Parameters::GetInstance(ParameterType::SM);
  double v1 = (*sm1)("GAUGE", LhaID(1)).real();
  //TODO : redo tests
  HyperisoConfig cfg2 = cfg;
  // cfg2.flags[ExternalFlag::USE_MARTY] = true;
  mm->switch_lha(lha2.string(), cfg2);

  sm1->CleanupInstance(ParameterType::SM);
  auto sm2 = Parameters::GetInstance(ParameterType::SM);

  double v2 = (*sm2)("GAUGE", LhaID(1)).real();
  std::cout << v1 << " " << v2 << std::endl;
  assert(v1 != v2);

  auto blist = sm2->get_blocks_list();
  assert(blist.find("MASS")  != blist.end());
  assert(blist.find("GAUGE") != blist.end());
  assert(blist.find("VCKM")  != blist.end());

  const auto& types = mm->getMemoryCache().parameter_types;
  assert(hasType(types, ParameterType::WILSON));

  std::cout << " INTEGRATION OK\n";
  return 0;
}

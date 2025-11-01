#include <cassert>
#include <iostream>
#include "Parameters.h"
#include "MemoryManager.h"
#include "Config.h"
#include "Paths.h"
#include "General.h"
#include "dummies_mm.hpp"

static bool contains(const std::unordered_set<BlockName>& s, const BlockName& b){
  return s.find(b) != s.end();
}

int main() {
  std::cout << "== Parameters UNIT ==\n";


    const fs::path lha_rel = "unit_params.lha";
    fs::path sandbox = prepare_assets_for_mm(lha_rel);

    auto ba   = std::make_shared<DummyBA>();
  auto cp   = std::make_shared<DummyCorr<ParamId>>();
  auto co   = std::make_shared<DummyCorr<ObservableId>>();
  auto spec = std::make_shared<DummySpectrum>();
  auto tpp = std::make_shared<TestPathsProvider>(sandbox);
  (void)MemoryManager::Create(ba, cp, co, spec, tpp);

  Config cfg;
  cfg.model = Model::SM;
  cfg.flags[ExternalFlag::HAS_WILSON_INPUT] = false;
  // cfg.flags[ExternalFlag::USE_MARTY] = false;

  MemoryManager::GetInstance()->init(lha_rel.string(), cfg);

  auto sm = Parameters::GetInstance(ParameterType::SM);

  assert(sm->exist("MASS", LhaID(25)));
  double v = (*sm)("MASS", LhaID(25));
  assert(v > 0.0);

  sm->setBlockValue("MASS", LhaID(25), v + 1.0);
  assert(std::abs((*sm)("MASS", LhaID(25)) - (v + 1.0)) < 1e-12);

  sm->shiftParameter(ParamId{ParameterType::SM, "MASS", LhaID(25)}, 2.0);
  assert(std::abs((*sm)("MASS", LhaID(25)) - (v + 3.0)) < 1e-12);

  sm->freeze_param("MASS", LhaID(25));
  sm->unfreeze_param("MASS", LhaID(25));
  sm->freeze_block("MASS");
  sm->unfreeze_block("MASS");

  double sc = sm->get_block_scale("MASS");
  assert(std::abs(sc - 91.1876) < 1e-6);

  auto blist = sm->get_blocks_list();
  assert(contains(blist, "MASS"));
  assert(contains(blist, "GAUGE"));
  assert(contains(blist, "VCKM"));

  std::cout << " UNIT OK\n";
  return 0;
}

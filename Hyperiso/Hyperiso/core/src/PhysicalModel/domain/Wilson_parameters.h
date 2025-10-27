#pragma once
#include <array>
#include "IBlockComposer.h"
#include "BWilsonRunningParameters.h"
#include "ParameterProxy.h"
#include "Include.h"
#include "IWilsonParameters.h"
#include "BWilsonRunningParameters.h"
#include "MesonMixingRunningParameters.h"

using BRP = BWilsonRunningParameters;
using MMRP = MesonMixingRunningParameters;

class WilsonParameterHelper : public IWilsonParameterHelper {
public:
    WilsonParameterHelper(std::shared_ptr<IBlockComposer> ibc) : IWilsonParameterHelper(ibc) {}
    void init(int gen, WGroupId grp) override;
    void cleanup() override;
    
private:
    void init_scale_independent_block(int gen) override;
    void init_matching_block() override;
    void init_running_block(WGroupId grp) override;

    void init_running_parameter_blocks_B();
    void init_running_parameter_blocks_MM();

};


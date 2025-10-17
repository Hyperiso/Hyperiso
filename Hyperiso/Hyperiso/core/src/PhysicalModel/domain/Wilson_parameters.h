#pragma once
#include <array>
#include "IBlockComposer.h"
#include "BWilsonRunningParameters.h"
#include "ParameterProxy.h"
#include "Include.h"
#include "IWilsonParameters.h"
#include "BWilsonRunningParameters.h"

using BRP = BWilsonRunningParameters;

class WilsonParameterHelper : public IWilsonParameterHelper {
public:
    WilsonParameterHelper(std::shared_ptr<IBlockComposer> ibc) : IWilsonParameterHelper(ibc) {}
    void init(int gen) override;
    void cleanup() override;
    
private:
    void init_scale_independent_block(int gen) override;
    void init_matching_block() override;
    void init_running_block() override;

    void init_running_parameter_blocks();

};


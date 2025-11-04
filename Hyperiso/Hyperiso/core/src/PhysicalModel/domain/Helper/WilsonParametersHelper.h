#ifndef WILSON_PARAMETER_HELPER_H
#define WILSON_PARAMETER_HELPER_H

#include <array>
#include "IBlockComposer.h"
#include "Include.h"
#include "IWilsonParameters.h"
#include "BWilsonRunningParameters.h"
#include "MesonMixingRunningParameters.h"
#include "QCDHelper.h"

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

#endif
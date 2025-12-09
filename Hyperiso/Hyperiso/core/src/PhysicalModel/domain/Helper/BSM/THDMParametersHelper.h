#ifndef THDM_PARAMETERS_HELPER_H
#define THDM_PARAMETERS_HELPER_H

#include <algorithm>
#include <array>
#include <functional>
#include "WilsonParamComposer.h"
#include "IWilsonParameters.h"
#include "Logger.h"
#include "SourcesView.h"
#include <iostream>

class thdm_parameters : public IWilsonParameterHelper {
public:

    thdm_parameters(std::shared_ptr<IBlockComposer> iblock_c) : IWilsonParameterHelper(iblock_c) {}
    void init(int gen, WGroupId grp) override;
    void cleanup() override {} //TODO
protected:
    void init_scale_independent_block(int gen) override;
    void init_matching_block() override;
    void init_running_block(WGroupId grp) override {}

};

#endif
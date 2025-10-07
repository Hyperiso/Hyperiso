#ifndef THDM_PARAMETERS_H
#define THDM_PARAMETERS_H

#include <algorithm>
#include <array>
#include <functional>
#include "WilsonParamComposer.h"
#include "IWilsonParameters.h"
#include "Logger.h"
#include <iostream>

class thdm_parameters : public IWilsonParameterHelper {
public:

    thdm_parameters(std::shared_ptr<IBlockComposer> iblock_c) : IWilsonParameterHelper(iblock_c) {}
    void init(int gen) override;
    void cleanup() override {} //TODO
protected:
    void init_scale_independent_block(int gen) override;
    void init_matching_block() override;
    void init_running_block() override {}
    // inline bool is_initialized() {return thdm_parameters::initialized;}
    // static inline WilsonParamComposer composer = WilsonParamComposer();
    // // static inline double current_mu_W{-1};
    // static inline bool initialized {false};

};

#endif
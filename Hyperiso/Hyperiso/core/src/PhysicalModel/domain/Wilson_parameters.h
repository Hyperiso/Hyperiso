#pragma once
#include <array>
#include "WilsonParamComposer.h"
#include "BWilsonRunningParameters.h"
#include "ParameterProxy.h"
#include "Include.h"

class WilsonParameterHelper {
public:
    static void init(int gen);
    static void set_gen(int gen);

private:
    static void init_scale_independent_block(int gen);
    static void init_matching_block();
    static void init_running_block();

    static inline WilsonParamComposer composer = WilsonParamComposer();
    static inline double current_mu_W{-1};
    static inline double current_mu_h{-1};
    static inline bool initialized {false};
};


#include <algorithm>
#include <array>
#include <functional>
#include "WilsonParamComposer.h"
#include "Logger.h"
#include <iostream>

class thdm_parameters {
public:

    static void init(double mu_W);
    static void init_scale_independant_block();
    static void init_matching_block(double mu_W);

    static inline WilsonParamComposer composer = WilsonParamComposer();
    // static inline double current_mu_W{-1};
    static inline bool initialized {false};

};
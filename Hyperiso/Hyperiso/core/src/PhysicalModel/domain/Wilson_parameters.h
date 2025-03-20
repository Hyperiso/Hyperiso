#pragma once
#include <array>
#include "WilsonParamComposer.h"
#include "BWilsonRunningParameters.h"

class WilsonParameters {
public:
    WilsonParameters(double mu_W, double mu_h, int gen);
    void set_mu_W(double mu_W);
    void set_mu(double mu);
    void set_gen(int gen);

private:
    WilsonParamComposer composer;

    void init_scale_independent_block(int gen);
    void init_matching_block(double mu_W);
    void init_running_block(double mu_W, double mu_h);
};


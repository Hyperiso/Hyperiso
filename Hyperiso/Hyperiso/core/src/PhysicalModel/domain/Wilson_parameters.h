#pragma once
#include <array>
#include "WilsonParamComposer.h"
#include "BWilsonRunningParameters.h"
#include "ModelParamAdapter.h"

struct WilsonRunningMatrix {
    std::array<std::array<double, BWilsonRunningParameters::array_size>, BWilsonRunningParameters::array_size> U0 {};
    std::array<std::array<double, BWilsonRunningParameters::array_size>, BWilsonRunningParameters::array_size> U1 {};
    std::array<std::array<double, BWilsonRunningParameters::array_size>, BWilsonRunningParameters::array_size> U2 {};
    std::array<std::array<double, BWilsonRunningParameters::array_size>, BWilsonRunningParameters::array_size> V0 {};
    std::array<std::array<double, BWilsonRunningParameters::array_size>, BWilsonRunningParameters::array_size> V1 {};

    std::array<double, BWilsonRunningParameters::array_size> etaMuPowers = {};
	std::array<double, BWilsonRunningParameters::array_size> etaMuPowers2 = {};
};

class WilsonParameters {
public:
    WilsonParameters(double mu_W, double mu_h, int gen);
    void set_mu_W(double mu_W);
    void set_mu(double mu);
    void set_gen(int gen);

    WilsonRunningMatrix get_matrix() {return this->w_run;}

private:
    WilsonParamComposer composer;

    void init_scale_independent_block(int gen);
    void init_matching_block(double mu_W);
    void init_running_block(double mu_W, double mu_h);

    WilsonRunningMatrix w_run;
};


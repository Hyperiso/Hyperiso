#pragma once
#include <array>
#include "WilsonParamComposer.h"
#include "BWilsonRunningParameters.h"
#include "ModelParamAdapter.h"
#include "Include.h"

struct WilsonRunningMatrix {
    std::array<std::array<double, BWilsonRunningParameters::array_size>, BWilsonRunningParameters::array_size> U0 {};
    std::array<std::array<double, BWilsonRunningParameters::array_size>, BWilsonRunningParameters::array_size> U1 {};
    std::array<std::array<double, BWilsonRunningParameters::array_size>, BWilsonRunningParameters::array_size> U2 {};
    std::array<std::array<double, BWilsonRunningParameters::array_size>, BWilsonRunningParameters::array_size> V0 {};
    std::array<std::array<double, BWilsonRunningParameters::array_size>, BWilsonRunningParameters::array_size> V1 {};

    std::array<double, BWilsonRunningParameters::array_size> etaMuPowers = {};
	std::array<double, BWilsonRunningParameters::array_size> etaMuPowers2 = {};
};

class BWilsonRunningHelper {
public:
    WilsonRunningMatrix get_matrix() {return this->w_run;}
    void update();

private:
    WilsonRunningMatrix w_run;
};

class WilsonParameterHelper {
public:
    static void init(double mu_W, double mu_h, int gen);
    static void set_mu_W(double mu_W);
    static void set_mu(double mu);
    static void set_gen(int gen);

private:
    static void init_scale_independent_block(int gen);
    static void init_matching_block(double mu_W);
    static void init_running_block(double mu_W, double mu_h);

    static inline WilsonParamComposer composer = WilsonParamComposer();
    static inline double current_mu_W{-1};
    static inline double current_mu_h{-1};
    static inline bool initialized {false};
};


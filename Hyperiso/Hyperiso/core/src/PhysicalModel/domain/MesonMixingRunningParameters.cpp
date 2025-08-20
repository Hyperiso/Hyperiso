#include "MesonMixingRunningParameters.h"

const std::array<std::array<double, MesonMixingRunningParameters::n_coefs>, MesonMixingRunningParameters::n_coefs> MesonMixingRunningParameters::SUSY_to_BMU_superiso {{
    {1.000, 0.000,  0.000, 0.000,  0.000, 0.000, 0.000,  0.000},
    {0.000, 0.000,  0.000, 0.000, -0.500, 0.000, 0.000,  0.000},
    {0.000, 0.000,  0.000, 1.000,  0.000, 0.000, 0.000,  0.000},
    {0.000, 1.000, -0.500, 0.000,  0.000, 0.000, 0.000,  0.000},
    {0.000, 0.000,  0.125, 0.000,  0.000, 0.000, 0.000,  0.000},
    {0.000, 0.000,  0.000, 0.000,  0.000, 1.000, 0.000,  0.000},
    {0.000, 0.000,  0.000, 0.000,  0.000, 0.000, 1.000, -0.500},
    {0.000, 0.000,  0.000, 0.000,  0.000, 0.000, 0.000,  0.125},
}};

const std::array<std::array<double, MesonMixingRunningParameters::n_coefs>, MesonMixingRunningParameters::n_coefs> MesonMixingRunningParameters::SUSY_to_BMU_fierz {{
    {1.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0},
    {0.0, 0.0,  0.0, 0.0, 2.5, 0.0, 0.0,  0.0},
    {0.0, 0.0,  0.0, 1.0, 0.0, 0.0, 0.0,  0.0},
    {0.0, 1.0,  2.5, 0.0, 0.0, 0.0, 0.0,  0.0},
    {0.0, 0.0, -1.5, 0.0, 0.0, 0.0, 0.0,  0.0},
    {0.0, 0.0,  0.0, 0.0, 0.0, 1.0, 0.0,  0.0},
    {0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 1.0,  2.5},
    {0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, -1.5},
}};

const std::array<std::array<std::array<double, 2>, 2>, 2> MesonMixingRunningParameters::a0_LR_5 {{
    {{{1.0000,  0.0000}, {0.0000, 0.0000}}},
    {{{0.6667, -0.6667}, {0.0000, 1.0000}}},
}};

const std::array<std::array<std::array<double, 2>, 2>, 2> MesonMixingRunningParameters::a1_LR_5 {{
    {{{ -2.0994,  0.9250}, {0.0000, -1.3875}}},
    {{{-11.7329, -5.3048}, {0.0000,  7.9572}}},
}};

const std::array<std::array<std::array<double, 2>, 2>, 2> MesonMixingRunningParameters::b_LR_5 {{
    {{{1.1744,  0.0000}, {1.3875,  0.0000}}},
    {{{0.7829, 16.2548}, {0.9250, -8.8822}}},
}};

const std::array<std::array<std::array<double, 2>, 2>, 2> MesonMixingRunningParameters::a0_S_5 {{
    {{{ 1.0153, -0.0153}, { 1.9325, -1.9325}}},
    {{{-0.0081,  0.0081}, {-0.0153,  1.0153}}},
}};

const std::array<std::array<std::array<double, 2>, 2>, 2> MesonMixingRunningParameters::a1_S_5 {{
    {{{4.8177,  0.3371}, {9.1696, 42.5021}}},
    {{{0.0531, -0.0566}, {0.1011, -7.1314}}},
}};

const std::array<std::array<std::array<double, 2>, 2>, 2> MesonMixingRunningParameters::b_S_5 {{
    {{{-5.2272,  0.0724}, {-38.8778, -12.7939}}},
    {{{ 0.0415,  0.3083}, { -0.0380,   6.7219}}},
}};

const std::array<std::array<double, 2>, 2> a0_LR_4_11 {{{1.0000,  0.0000}, {0.0000,  0.0000}}};
const std::array<std::array<double, 2>, 2> a0_LR_4_12 {{{0.0000,  0.0000}, {0.0000,  0.0000}}};
const std::array<std::array<double, 2>, 2> a0_LR_4_21 {{{0.6667,  0.0000}, {0.0000, -0.6667}}};
const std::array<std::array<double, 2>, 2> a0_LR_4_22 {{{0.0000,  0.0000}, {0.0000,  1.0000}}};
const std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> MesonMixingRunningParameters::a0_LR_4 {{{a0_LR_4_11, a0_LR_4_12}, {a0_LR_4_21, a0_LR_4_22}}};

const std::array<std::array<double, 2>, 2> a1_LR_4_11 {{{ -2.0241,  0.0000}, {0.0000,  0.9279}}};
const std::array<std::array<double, 2>, 2> a1_LR_4_12 {{{  0.0000,  0.0000}, {0.0000, -1.3918}}};
const std::array<std::array<double, 2>, 2> a1_LR_4_21 {{{-16.6828,  0.0000}, {0.0000, -4.4701}}};
const std::array<std::array<double, 2>, 2> a1_LR_4_22 {{{  0.0000,  0.0000}, {0.0000,  6.7052}}};
const std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> MesonMixingRunningParameters::a1_LR_4 {{{a1_LR_4_11, a1_LR_4_12}, {a1_LR_4_21, a1_LR_4_22}}};

const std::array<std::array<double, 2>, 2> b_LR_4_11 {{{-0.0753, -0.0029}, {0.0000,  0.0000}}};
const std::array<std::array<double, 2>, 2> b_LR_4_12 {{{ 0.0000,  0.0043}, {0.0000,  0.0000}}};
const std::array<std::array<double, 2>, 2> b_LR_4_21 {{{-0.0502, -0.0019}, {0.0000,  0.0029}}};
const std::array<std::array<double, 2>, 2> b_LR_4_22 {{{ 5.0000, -0.8327}, {0.0000,  1.2491}}};
const std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> MesonMixingRunningParameters::b_LR_4 {{{b_LR_4_11, b_LR_4_12}, {b_LR_4_21, b_LR_4_22}}};

const std::array<std::array<double, 2>, 2> c_LR_4_11 {{{1.1744,  0.0000}, {0.0000,  0.0000}}};
const std::array<std::array<double, 2>, 2> c_LR_4_12 {{{1.3875,  0.0000}, {0.0000,  0.0000}}};
const std::array<std::array<double, 2>, 2> c_LR_4_21 {{{0.7829,  0.0000}, {0.0000, 16.2548}}};
const std::array<std::array<double, 2>, 2> c_LR_4_22 {{{0.9250,  0.0000}, {0.0000, -8.8822}}};
const std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> MesonMixingRunningParameters::c_LR_4 {{{c_LR_4_11, c_LR_4_12}, {c_LR_4_21, c_LR_4_22}}};

const std::array<std::array<double, 2>, 2> a0_S_4_11 {{{ 1.0153,  0.0000}, {0.0000, -0.0153}}};
const std::array<std::array<double, 2>, 2> a0_S_4_12 {{{ 1.9325,  0.0000}, {0.0000, -1.9325}}};
const std::array<std::array<double, 2>, 2> a0_S_4_21 {{{-0.0081,  0.0000}, {0.0000,  0.0081}}};
const std::array<std::array<double, 2>, 2> a0_S_4_22 {{{-0.0153,  0.0000}, {0.0000,  1.0153}}};
const std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> MesonMixingRunningParameters::a0_S_4 {{{a0_S_4_11, a0_S_4_12}, {a0_S_4_21, a0_S_4_22}}};

const std::array<std::array<double, 2>, 2> a1_S_4_11 {{{ 4.2458,  0.0000}, {0.0000,  0.3640}}};
const std::array<std::array<double, 2>, 2> a1_S_4_12 {{{ 8.0810,  0.0000}, {0.0000, 45.9008}}};
const std::array<std::array<double, 2>, 2> a1_S_4_21 {{{ 0.0587,  0.0000}, {0.0000, -0.0534}}};
const std::array<std::array<double, 2>, 2> a1_S_4_22 {{{ 0.1117,  0.0000}, {0.0000, -6.7398}}};
const std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> MesonMixingRunningParameters::a1_S_4 {{{a1_S_4_11, a1_S_4_12}, {a1_S_4_21, a1_S_4_22}}};

const std::array<std::array<double, 2>, 2> b_S_4_11 {{{ 0.5700, -0.0334}, { 0.0020,  0.0064}}};
const std::array<std::array<double, 2>, 2> b_S_4_12 {{{ 1.0848, -4.2075}, { 0.0038,  0.8087}}};
const std::array<std::array<double, 2>, 2> b_S_4_21 {{{-0.0045,  0.0003}, {-0.0011, -0.0034}}};
const std::array<std::array<double, 2>, 2> b_S_4_22 {{{-0.0086,  0.0334}, {-0.0020, -0.4249}}};
const std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> MesonMixingRunningParameters::b_S_4 {{{b_S_4_11, b_S_4_12}, {b_S_4_21, b_S_4_22}}};

const std::array<std::array<double, 2>, 2> c_S_4_11 {{{ -5.2272, 0.0000}, {0.0000,   0.0724}}};
const std::array<std::array<double, 2>, 2> c_S_4_12 {{{-38.8778, 0.0000}, {0.0000, -12.7939}}};
const std::array<std::array<double, 2>, 2> c_S_4_21 {{{  0.0415, 0.0000}, {0.0000,  -0.0380}}};
const std::array<std::array<double, 2>, 2> c_S_4_22 {{{  0.3083, 0.0000}, {0.0000,   6.7219}}};
const std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> MesonMixingRunningParameters::c_S_4 {{{c_S_4_11, c_S_4_12}, {c_S_4_21, c_S_4_22}}};
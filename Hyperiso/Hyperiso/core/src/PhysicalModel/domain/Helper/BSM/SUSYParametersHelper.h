#ifndef SUSY_PARAMETERS_HELPER_H
#define SUSY_PARAMETERS_HELPER_H

#include <algorithm>
#include <array>
#include <functional>
// #include "epsilon_calculator.h"
#include "Math_SUSY.h"
#include "Logger.h"
#include "ParameterProxy.h"
#include "WilsonParamComposer.h"
#include "IWilsonParameters.h"

constexpr double Pi = 3.14159265358979323846;

constexpr int N_UL_UR = 7, M_UL_UR = 4, N_NL_NR = 4, M_NL_NR = 4, N_Gamma_U = 7, M_Gamma_U = 7, N_X = 3, N_Mch = 3, N_MsqU = 6, N_MsqD = 6, N_Msn = 4;

using Array2D_7x4 = std::array<std::array<double, M_UL_UR>, N_UL_UR>;
using Array2D_4x4 = std::array<std::array<double, M_NL_NR>, N_NL_NR>;
using Array2D_7x7 = std::array<std::array<double, M_Gamma_U>, N_Gamma_U>;
using Array3D_3x7x4 = std::array<Array2D_7x4, N_X>;
using Array1D_4 = std::array<double, 4>;
using Array1D_3 = std::array<double, N_Mch>;
using Array1D_7 = std::array<double, N_MsqU>;
using Array2D_4x4_I = std::array<std::array<complex_t, M_NL_NR>, N_NL_NR>;


class susy_parameters : public IWilsonParameterHelper {
public:

	susy_parameters(std::shared_ptr<IBlockComposer> iblock_c) : IWilsonParameterHelper(iblock_c) {}

	void init(int gen, WGroupId grp) override;
	void cleanup() override {}

protected:
	void init_epsilon_block();
    void init_scale_independent_block(int gen) override;
    void init_matching_block() override;
	void init_running_block(WGroupId grp) override {};
};

#endif
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

// struct WilsonSusyMatrix {
// 	Array2D_7x4 Gamma_UL {};
// 	Array2D_7x4 Gamma_UR {};
// 	Array2D_4x4 Gamma_NL {};
// 	Array2D_4x4 Gamma_NR {};
// 	Array2D_7x7 Gamma_U{};
// 	Array2D_7x7 I_LR{};
// 	Array2D_7x7 P_U{};
// 	Array3D_3x7x4 X_UL{};
// 	Array3D_3x7x4 X_UR{};
// 	Array3D_3x7x4 X_NL{};
// 	Array3D_3x7x4 X_NR{};
// 	std::array<std::array<std::array<std::array<double, 4>, 4>, 3>, 7> G_aimn;

// };

class susy_parameters : public IWilsonParameterHelper {
	// bool is_PrimeCQG = false;

    // explicit susy_parameters();
    // susy_parameters(const susy_parameters&) = delete;
    // void operator=(const susy_parameters&) = delete;

public:

	susy_parameters(std::shared_ptr<IBlockComposer> iblock_c) : IWilsonParameterHelper(iblock_c) {}
	// void update();

	void init(int gen, WGroupId grp) override;
	void cleanup() override {}

protected:
	void init_epsilon_block();
    void init_scale_independent_block(int gen) override;
    void init_matching_block() override;
	void init_running_block(WGroupId grp) override {};
	// void reset_PrimeCQG(double Qmatch);
	// void reset_G();
	
	// Array2D_4x4  VCKM = {{{0.,0.,0.,0.},
	// 				{param->Vud, param->Vus, -(param->Vts * param->Vtb + param->Vcs * param->Vcb) / param->Vus, 0.0},
	// 				{param->Vcd, param->Vcs, param->Vcb, 0.0},
	// 				{param->Vtd, param->Vts, param->Vtb, 0.0}
	// 				}}; 
	// Array2D_4x4_I  VCKM = {{{(*sm)("CKM",11), (*sm)("CKM",12), -((*sm)("CKM",32) * (*sm)("CKM",33) + (*sm)("sm",22) * (*sm)("CKM",23)) / (*sm)("CKM",12)},
	// 				{(*sm)("CKM",21), (*sm)("CKM",22), (*sm)("CKM",23)},
	// 				{(*sm)("CKM",31), (*sm)("CKM",32),(*sm)("CKM",33)}
	// 				}}; 

	// complex_t c11 = (*sm)("RECKM", 00) + (*sm)("IMCKM", 00) * complex_t(0, 1);
	// complex_t c12 = (*sm)("RECKM", 01) + (*sm)("IMCKM", 01) * complex_t(0, 1);
	// complex_t c13 = (*sm)("RECKM", 02) + (*sm)("IMCKM", 02) * complex_t(0, 1);
	// complex_t c21 = (*sm)("RECKM", 10) + (*sm)("IMCKM", 10) * complex_t(0, 1);
	// complex_t c22 = (*sm)("RECKM", 11) + (*sm)("IMCKM", 11) * complex_t(0, 1);
	// complex_t c23 = (*sm)("RECKM", 12) + (*sm)("IMCKM", 12) * complex_t(0, 1);
	// complex_t c31 = (*sm)("RECKM", 20) + (*sm)("IMCKM", 20) * complex_t(0, 1);
	// complex_t c32 = (*sm)("RECKM", 21) + (*sm)("IMCKM", 21) * complex_t(0, 1);
	// complex_t c33 = (*sm)("RECKM", 22) + (*sm)("IMCKM", 22) * complex_t(0, 1);

	
	// complex_t complexTerm = -(c32 * c33 + c22 * c23) / c12;

	// Array2D_4x4_I VCKM = {{
    // {
    //     c11,
    //     c12,
    //     c13
    // },
    // {
    //     c21,
    //     c22,
    //     c23
    // },
    // {
    //     c31,
    //     c32,
    //     c33
    // }
// }};
	



    // double mass_H03 = (*susy)("MASS",36); // We have H_0^3 = A_0 from pdg numering scheme ?
    // double mass_A02 = (*susy)("HMIX",4);


    // double kappaFactor;
    // double B0c1 = 0.0, B0c2 = 0.0, B90c = 0.0, B100c = 0.0, C90c = 0.0, D90c = 0.0;
    // bool test;
    

	// WilsonSusyMatrix w_susy;
};

#endif
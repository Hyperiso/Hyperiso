#ifndef __BWILSONRUNNINGPARAMETERS_H__
#define __BWILSONRUNNINGPARAMETERS_H__

#include <array>

struct BWilsonRunningParameters {
    static constexpr int array_size {10};

    static const std::array<std::array<std::array<double, array_size>, array_size>, array_size> m00;
    static const std::array<std::array<std::array<double, array_size>, array_size>, array_size> m10;
    static const std::array<std::array<std::array<double, array_size>, array_size>, array_size> m11;
    static const std::array<std::array<std::array<double, array_size>, array_size>, array_size> m20;
    static const std::array<std::array<std::array<double, array_size>, array_size>, array_size> m21;
    static const std::array<std::array<std::array<double, array_size>, array_size>, array_size> m22;

    static const std::array<std::array<std::array<double, array_size>, array_size>, array_size> l00;
    static const std::array<std::array<std::array<double, array_size>, array_size>, array_size> l01;
    static const std::array<std::array<std::array<double, array_size>, array_size>, array_size> l10;
    static const std::array<std::array<std::array<double, array_size>, array_size>, array_size> l11;

    static constexpr std::array<double, array_size> ai = {14.0 / 23.0, 16.0 / 23.0, 6.0 / 23.0, -12.0 / 23.0, 0.408619, -0.422989, -0.899395, 0.145649, -1.0, -1.0};
    static constexpr std::array<double, array_size> ai2 = {6./23., -12./23., 0.4086, -0.4230, -0.8994, 0.1456, 16./23., 14./23., 11./23., 29./23.}; 

    static constexpr double a[8]={-2.,-1.,-0.899395,-0.521739,-0.422989,0.145649,0.260870,0.408619};
	static constexpr double b[8]={0.00354,0.01223,-0.00977,-0.01070,-0.00572,0.00022,0.01137,-0.00117};
	static constexpr double d_2a[8]={0.,0.,0.61602,0.44627,0.57472,0.08573,-0.48807,-0.24089};
	static constexpr double d_2b[8]={-1.18162,0.22940,0.06522,-0.04380,-0.02201,-0.00316,-0.03366,-0.00414};
	static constexpr double d_1[8]={0.01117,-0.03088,0.00411,0.00713,0.00478,0.00012,0.00379,-0.00023};
	static constexpr double d_4[8]={-0.00799,-0.03666,0.06300,0.,-0.01519,-0.00071,0.,-0.00344};
	static constexpr double e_1a[8]={0.,0.,-0.25941,-0.29751,-0.48014,0.04647,-0.16269,-0.04728};
	static constexpr double e_1b[8]={1.13374,0.09381,-0.03041,0.00781,0.01838,-0.00138,-0.02259,0.00121};
	static constexpr double e_4a[8]={0.,0.,-4.03683,0.,1.52565,-0.27461,0.,-0.70642};
	static constexpr double e_4b[8]={3.38669,-0.10885,0.16283,0.,0.06697,-0.01681,0.,0.00137,};
	static constexpr double e_1[8]={0.01117,-0.03088,0.00411,0.00713,0.00478,0.00012,0.00379,-0.00023};
	static constexpr double e_2[8]={0.00354,0.01223,-0.00977,-0.01070,-0.00572,0.00022,0.01137,-0.00117};
	static constexpr double e_3[8]={0.02179,-0.12336,0.07870,0.,0.01930,0.00873,0.,-0.00516};
	static constexpr double e_4[8]={-0.00799,-0.03666,0.06400,0.,-0.01519,-0.00071,0.,-0.00344};
	static constexpr double e_5[8]={0.19550,-0.93249,0.37858,0.,0.39909,0.05921,0.,-0.09989};
	static constexpr double e_6[8]={-0.17154,0.39616,0.01201,0.,-0.19423,0.00357,0.,-0.04597};

	static constexpr double y_std[8] = {0, 0, -1./3, -4./9, -20./3, -80./9, 1, 0};
	static constexpr double z_std[8] = {0, 0, 1, -1./6, 20, -10./3, 0, 1};

	static inline complex_t C7_eff_std(const std::array<complex_t, 10>& Ci) {
		complex_t C {0};
            for (size_t k = 0; k < 8; k++) {
                C += Ci[k] * BWilsonRunningParameters::y_std[k];
            }
            return C;
	}

	static inline complex_t C8_eff_std(const std::array<complex_t, 10>& Ci) {
		complex_t C {0};
            for (size_t k = 0; k < 8; k++) {
                C += Ci[k] * BWilsonRunningParameters::z_std[k];
            }
            return C;
	}

	static const std::array<std::array<double, 8>, 8> std_to_trad_LO;
	static const std::array<std::array<double, 8>, 8> std_to_trad_NLO;

	static constexpr double exp_prime_running[12] = {0, 0, 0, 0, 0, 0, 16./23, 14./23, 0, 0, -12./23, -12./23};
};

#endif // __BWILSONRUNNINGPARAMETERS_H__

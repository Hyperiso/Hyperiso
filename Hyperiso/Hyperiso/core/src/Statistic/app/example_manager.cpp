#include "StatisticManager.h"
#include "ObservableInterfaceAdapter2.h"
#include "StatCorrelationProxy.h"
#include "StatParameterProxy.h"
#include "ObservableInterface.h"
#include "StatParamSourcesProxy.h"
#include <cassert>

void print_vec(const std::vector<double>& vec) {
    std::cout << "[ ";
    for (size_t i = 0; i < vec.size(); i++) {
        std::cout << vec.at(i) << (i == vec.size() - 1 ? " " : ", ");
    }
    std::cout << "]" << std::endl;
}

int main(int argc, char** argv) {
    // RealValuedForm F_gauss_nll = [](const Vec& p) {
    //     // mean
    //     double mu0 = 1e-12;
    //     double mu1 = 1e3;

    //     // covariance
    //     double s0 = 1e-13;
    //     double s1 = 1e2;
    //     double rho = 0.3;

    //     double det = s0*s0*s1*s1*(1 - rho*rho);

    //     // inverse covariance
    //     double a =  1.0/(s0*s0*(1-rho*rho));
    //     double b = -rho/(s0*s1*(1-rho*rho));
    //     double c =  1.0/(s1*s1*(1-rho*rho));

    //     double d0 = p[0] - mu0;
    //     double d1 = p[1] - mu1;

    //     return 0.5 * (a*d0*d0 + 2*b*d0*d1 + c*d1*d1);
    // };

    // MinimizationContext ctx;
    // ctx.max_iter = 500;
    // ctx.step_size = 0.1;
    // ctx.tol = 1e-4;
    // MinimizationResult mr = minimize_NM(F_gauss_nll, {3e-12, 250}, {1e-13, 1e2}, ctx);
    // std::cout << mr.min << std::endl;
    // std::cout << mr.argmin[0] << ", " << mr.argmin[1] << std::endl;
    
    // Vec p_hat = {1e-12, 1e3};
    // auto H = hessian(F_gauss_nll, p_hat);
    // auto H_inv = inverse_hessian(F_gauss_nll, p_hat);

    // std::cout << H << std::endl;
    // std::cout << H_inv << std::endl;
    // std::cout << H * H_inv << std::endl;

    HyperisoMaster hyp;
    HyperisoConfig config_hyp;
    config_hyp.model = Model::SM;

    LOG_INFO("YO");

    hyp.init("lha/si_input.flha", config_hyp);

    StatisticConfig config;
    config.MC_draws = 100;
    config.obss = {
        {ObservableMapper::to_id(Observables::BR_BS_MUMU_UNTAG), QCDOrder::LO},
        {ObservableMapper::to_id(Observables::BR_BD_MUMU), QCDOrder::LO}
    };
    config.MLE_max_iter = 1000;
    config.MLE_tol = 1e-6;

    config.p_specs = {
        ParamId(ParameterType::SM, "VCKMIN", 1),
        ParamId(ParameterType::SM, "VCKMIN", 2),
        ParamId(ParameterType::SM, "VCKMIN", 3),
        ParamId(ParameterType::SM, "VCKMIN", 4)
    };

    std::shared_ptr<ObservableInterface> oi = std::make_shared<ObservableInterface>();

    StatisticManager stat(config, std::make_shared<ObservableInterfaceAdapterObs>(oi), std::make_shared<StatCorrelationProxy>(), std::make_shared<StatParameterProxy>(), std::make_shared<StatParamSourcesProxy>());
    LOG_INFO("YO1");
    stat.fill_cache();
    // LOG_INFO("YO2");
    auto start = std::chrono::steady_clock::now();
    // stat.compute_uncertainties();
    auto stop  = std::chrono::steady_clock::now();
    LOG_INFO("YO3");
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    std::cout << "Uncertainty estimation time : " << us << " µs\n";

    // auto start = std::chrono::steady_clock::now();
    auto fr = stat.compute_MLE();
    // auto stop  = std::chrono::steady_clock::now();
    // auto us = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    std::cout << "MLE fitting time : " << us << " µs\n";

    auto p_hat2 = fr.p_hat;
    auto p_hat_std = fr.p_hat_std;

    // for (auto [k, v] : p_hat) {
    //     std::cout << k << " = " << v << "(" << p_hat_std.at(k) << ")" << std::endl;
    // }

    std::cout << fr.p_correlations << std::endl;

    // auto paths = stat.confidence_contour(
    //     ParamId(ParameterType::SM, "VCKMIN", 3),
    //     ParamId(ParameterType::SM, "VCKMIN", 4),
    //     1,
    //     {{0.13, 0.19, 0.34, 0.38}}
    // );
    
    // for (auto &p : paths) {
    //     for (auto& xy: p) {
    //         std::cout << xy.first << " " << xy.second << std::endl;
    //     }
    //     std::cout << std::endl;
    // }

    return 0;
}

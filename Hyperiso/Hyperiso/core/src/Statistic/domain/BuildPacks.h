#pragma once
#include <vector>
#include "LinearAlgebra.h"
#include "StatParameterProxy.h"
#include "StatCorrelationProxy.h"
#include "ObservableInterface.h" #TODO

struct ExpPack { Vec mean; Matrix cov; };
struct NuisancePack { Vec mean; Matrix cov; };


inline ExpPack build_exp_from_proxies(const std::vector<Observables>& obs,
        CorrelationProvider::CorrelationType ctype,
        UncertaintyType u_type = UncertaintyType::COMBINED) {
    ExpPack P; const std::size_t n = obs.size();
    P.mean.resize(n); P.cov.assign(n, Vec(n, 0.0));
    StatCorrelationProxy cp;


    std::vector<double> sigma(n);
    for (std::size_t i=0;i<n;++i) {
        ObservableInterface oi;
        P.mean[i] = oi.get_exp_value(obs[i]);
        sigma[i] = oi.get_exp_uncertainty(obs[i], u_type);
    }
    for (std::size_t i=0;i<n;++i) for (std::size_t j=0;j<=i;++j) {
        double rho = cp(obs[i], obs[j], ctype);
        P.cov[i][j] = P.cov[j][i] = rho * sigma[i] * sigma[j];
    }
    return P;
}


inline NuisancePack build_eta_from_param_proxies(const std::vector<ParamId>& ids,
        CorrelationProvider::CorrelationType ctype,
        ParameterType ptype = ParameterType::SM) {
    const std::size_t n = ids.size();
    NuisancePack P; P.mean.resize(n); P.cov.assign(n, Vec(n, 0.0));
    StatParameterProxy pp(ptype);
    StatCorrelationProxy cp;
    std::vector<double> sigma(n);
    for (std::size_t i=0;i<n;++i) {
        P.mean[i] = pp(ids[i], DataType::VALUE);
        sigma[i] = pp(ids[i], DataType::STD_COMBINED);
    }
    for (std::size_t i=0;i<n;++i) for (std::size_t j=0;j<=i;++j) {
        double rho = cp(ids[i], ids[j], ctype);
        P.cov[i][j] = P.cov[j][i] = rho * sigma[i] * sigma[j];
    }
    return P;
}
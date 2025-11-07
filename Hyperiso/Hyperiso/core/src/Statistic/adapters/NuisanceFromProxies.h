#pragma once
#include <vector>
#include <stdexcept>
#include <unordered_map>
#include "StatParameterProxy.h"
#include "StatCorrelationProxy.h"
#include "RNGHelper.h"

// Helper that builds η̄ and Σ_η from proxies for a list of ParamId
struct NuisancePack {
Vector mean; Matrix cov;
};


inline NuisancePack build_nuisance_from_proxies(const std::vector<ParamId>& ids,
CorrelationProvider::CorrelationType ctype,
ParameterType ptype = ParameterType::SM) {
const std::size_t n = ids.size();
NuisancePack pack; pack.mean.resize(n); pack.cov.assign(n, Vector(n, 0.0));


StatParameterProxy pp(ptype);
StatCorrelationProxy cp;


for (std::size_t i=0;i<n;++i) pack.mean[i] = pp(ids[i], DataType::VALUE);


// Fill covariance from (σ_i, ρ_ij)
std::vector<double> sigma(n);
for (std::size_t i=0;i<n;++i) sigma[i] = pp(ids[i], DataType::STD_COMBINED); //USE combined here TODO


for (std::size_t i=0;i<n;++i) {
for (std::size_t j=0;j<=i;++j) {
const double rho = cp(ids[i], ids[j], ctype);
const double cij = rho * sigma[i] * sigma[j];
pack.cov[i][j] = pack.cov[j][i] = cij;
}
}
return pack;
}
#ifndef COVARIANCE_TRANSFORMER_H
#define COVARIANCE_TRANSFORMER_H

#include "IStatParameterProxy.h"
#include "General.h"

class CovarianceTransformer {
public:
    CovarianceTransformer(IStatParameterProxy<std::string, LhaID> corr_proxy) : corr_proxy(corr_proxy) { }


    std::vector<std::vector<double>> transform(std::vector<ParamId> ids);
    std::vector<std::vector<double>> transform(std::vector<ObservableId> ids);

private:
    IStatParameterProxy<std::string, LhaID> corr_proxy;
};

#endif
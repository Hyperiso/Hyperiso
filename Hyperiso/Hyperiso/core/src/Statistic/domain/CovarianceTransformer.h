#ifndef COVARIANCE_TRANSFORMER_H
#define COVARIANCE_TRANSFORMER_H

#include "IStatParameterProxy.h"
#include "General.h"
#include "CorrelationProvider.h" //For enum
#include "IStatCorrelationProxy.h"
#include "StatParameterProxy.h"

class CovarianceTransformer {
public:
    CovarianceTransformer(std::shared_ptr<IStatCorrelationProxy> corr_proxy, std::shared_ptr<IStatParameterProxy> par_proxy) : corr_proxy(corr_proxy) { }


    std::vector<std::vector<double>> transform(std::vector<ParamId> ids);
    std::vector<std::vector<double>> transform(std::vector<ObservableId> ids);

private:
    std::shared_ptr<IStatCorrelationProxy> corr_proxy;
};

#endif
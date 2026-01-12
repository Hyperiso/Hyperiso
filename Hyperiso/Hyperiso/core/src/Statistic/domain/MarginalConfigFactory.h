#ifndef __MARGINALCONFIGFACTORY_H__
#define __MARGINALCONFIGFACTORY_H__

#include "Include.h"
#include "MarginalType.h"
#include "StatParameterProxy.h"
#include "GaussianMarginal.h"
#include "FlatMarginal.h"
#include "SplitGaussianMarginal.h"
#include "LikelihoodMarginal.h"

using Config = std::variant<FlatMarginalCfg, GaussianMarginalCfg, SplitGaussianMarginalCfg, LikelihoodMarginalCfg>;

class MarginalConfigFactory {
public:
    static Config create(ParamId pid, MarginalType marginal);

private:
    static inline const StatParameterProxy p {};
};

#endif // __MARGINALCONFIGFACTORY_H__

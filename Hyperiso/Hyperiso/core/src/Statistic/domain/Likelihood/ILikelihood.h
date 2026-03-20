#ifndef __ILIKELIHOOD_H__
#define __ILIKELIHOOD_H__

#include "Math.h"

class ILikelihood {
public:
    virtual ~ILikelihood() = default;
    virtual double nll(const std::vector<double>& theta) const = 0;
    virtual std::vector<fit_app::ParameterDefinition> get_param_defs() const = 0;
    virtual std::size_t dim() const = 0;
};

#endif // __ILIKELIHOOD_H__

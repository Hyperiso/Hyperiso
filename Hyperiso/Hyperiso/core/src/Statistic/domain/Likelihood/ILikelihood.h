#ifndef __ILIKELIHOOD_H__
#define __ILIKELIHOOD_H__

#include "Math.h"

class ILikelihood {
public:
    virtual ~ILikelihood() = default;
    virtual double nll(const Vector& theta) const = 0;
    virtual std::vector<fit_app::ParameterDefinition> get_param_defs() const = 0;

    std::size_t dim() const { return this->dim_; }

protected:
    std::size_t dim_;
};

#endif // __ILIKELIHOOD_H__

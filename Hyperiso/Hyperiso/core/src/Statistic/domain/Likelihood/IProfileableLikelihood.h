#ifndef IPROFILEABLE_LIKELIHOOD_H
#define IPROFILEABLE_LIKELIHOOD_H

#include "ILikelihood.h"

class IProfileableLikelihood : public ILikelihood {
public:
    virtual std::size_t p_dimension() const = 0;
    virtual std::size_t eta_dimension() const = 0;

    virtual std::vector<double> central_p() const = 0;
    virtual std::vector<double> central_eta() const = 0;

    virtual std::vector<double> predict(
        const std::vector<double>& p,
        const std::vector<double>& eta
    ) const = 0;

    virtual std::vector<double> residuals(
        const std::vector<double>& p,
        const std::vector<double>& eta
    ) const = 0;

    virtual double nll_from_split(
        const std::vector<double>& p,
        const std::vector<double>& eta
    ) const = 0;

    virtual RealMatrix observable_curvature(
        const std::vector<double>& residuals
    ) const = 0;

    virtual RealMatrix nuisance_curvature(
        const std::vector<double>& eta
    ) const = 0;
};

#endif
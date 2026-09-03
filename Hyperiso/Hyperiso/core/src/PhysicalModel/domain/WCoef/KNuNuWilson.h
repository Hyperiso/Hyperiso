#ifndef KNUNU_WILSON_H
#define KNUNU_WILSON_H

#include "Wilson.h"

/**
 * @brief Flavour-resolved diagonal WET coefficients for s -> d nu_i anti-nu_i.
 *
 * K->pi nu anti-nu is evaluated exclusively with this flavour-resolved basis.
 * The native SM matching functions are shared with the historical CK_L
 * implementation to avoid duplicating the Inami-Lim calculation, but the decay
 * no longer falls back to a flavour-universal coefficient.
 */
class KNuNuWilson : public WilsonCoefficient {
public:
    KNuNuWilson(WCoef coefficient, ContributionType contribution);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<KNuNuWilson>(*this);
    }

    static bool is_left(WCoef coefficient);
    static bool is_sm_diagonal_left(WCoef coefficient);

private:
    WCoef coefficient_;
};

#endif

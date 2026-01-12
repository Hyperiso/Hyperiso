#ifndef IDISTRIBUTION_H
#define IDISTRIBUTION_H

#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include "RNGHelper.h"

struct IMarginalDistribution {
    virtual ~IMarginalDistribution() = default;

    virtual Vector rvs(std::size_t n) = 0;
    virtual double logpdf(double x) = 0;
    virtual double cdf(double x) = 0;
    virtual double ppf(double x) = 0;
    virtual double mean() = 0;
    virtual double std() = 0;
};

#endif
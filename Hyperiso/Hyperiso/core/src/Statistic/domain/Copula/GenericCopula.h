#ifndef __GENERICCOPULA_H__
#define __GENERICCOPULA_H__

#include "ICopula.h"
#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

constexpr double CLIP_U {1e-15};

class GenericCopula: public ICopula {
public:
    explicit GenericCopula(unsigned int seed = std::random_device{}());

protected:
    const gsl_rng_type* rng_tp {gsl_rng_mt19937};
    gsl_rng* eng_;
};

#endif // __GENERICCOPULA_H__

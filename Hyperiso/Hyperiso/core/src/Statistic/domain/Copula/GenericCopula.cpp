#include "GenericCopula.h"

GenericCopula::GenericCopula(unsigned int seed)
    : eng_(gsl_rng_alloc(rng_tp), &gsl_rng_free)
{
    gsl_rng_set(eng_.get(), seed);
}
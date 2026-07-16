#ifndef GENERICCOPULA_H
#define GENERICCOPULA_H

#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <memory>

#include "ICopula.h"

using gsl_rng_sptr = std::unique_ptr<gsl_rng, decltype(&gsl_rng_free)>;

/**
 * @file GenericCopula.h
 * @brief Base class providing RNG infrastructure for concrete copulas.
 *
 * This class does not implement the @ref ICopula interface fully by itself,
 * but centralizes the random-number generator setup shared by concrete copula
 * implementations.
 *
 * It uses a GSL RNG engine (currently Mersenne Twister) seeded at construction.
 *
 * @see ICopula
 * @see GaussianCopula
 * @see StudentTCopula
 */

/// Small clipping threshold used to avoid exact 0 or 1 uniforms.
constexpr double CLIP_U {1e-15};

/**
 * @class GenericCopula
 * @brief Partial base implementation of a copula with RNG support.
 *
 * This class mainly exists to factorize:
 * - GSL RNG type selection,
 * - generator allocation and initialization.
 *
 * Derived classes are responsible for:
 * - specifying the dependence structure,
 * - implementing sampling and log-density evaluation.
 */
class GenericCopula: public ICopula {
public:
    /**
     * @brief Constructs the copula base with a seeded random engine.
     *
     * @param seed Seed used to initialize the underlying GSL RNG.
     */
    explicit GenericCopula(unsigned int seed = std::random_device{}());
    
protected:
    const gsl_rng_type* rng_tp {gsl_rng_mt19937};   /// GSL RNG type used by the copula.
    gsl_rng_sptr eng_{nullptr, &gsl_rng_free};      /// unique ptr of GSL random engine used by derived classes.
};

#endif // GENERICCOPULA_H

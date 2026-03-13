#ifndef ICOPULA_H
#define ICOPULA_H

#include "RNGHelper.h"

/**
 * @file ICopula.h
 * @brief Abstract interface for multivariate copulas.
 *
 * A copula models the dependence structure between random variables
 * independently of their marginal distributions.
 *
 * Implementations of this interface are expected to work on the unit cube:
 * - sampling returns vectors in \f$[0,1]^d\f$,
 * - density evaluation is performed on pseudo-observations / uniforms.
 *
 * Typical workflow:
 * @code
 *   std::unique_ptr<ICopula> cop = ...;
 *   Vector u = cop->sample_u();
 *   double logc = cop->log_density(u);
 * @endcode
 *
 * @see GaussianCopula
 * @see StudentTCopula
 */
class ICopula {
public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~ICopula() = default;

    /**
     * @brief Draws multiple samples from the copula.
     *
     * Each sample is a vector of uniform variates on \f$[0,1]\f$ whose
     * dependence is governed by the copula.
     *
     * @param n Number of samples to generate.
     * @return A vector of samples, each of dimension equal to the copula dimension.
     */
    virtual std::vector<Vector> sample_u(std::size_t n) = 0;

    /**
     * @brief Draws a single sample from the copula.
     *
     * @return One dependent uniform sample in \f$[0,1]^d\f$.
     */
    virtual Vector sample_u() = 0;

    /**
     * @brief Evaluates the log-density of the copula at a point in the unit cube.
     *
     * @param u Point in \f$[0,1]^d\f$.
     * @return \f$\log c(u)\f$, where \f$c\f$ is the copula density.
     */
    virtual double log_density(Vector u) = 0;
};

#endif // ICOPULA_H

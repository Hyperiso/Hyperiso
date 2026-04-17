#ifndef IDISTRIBUTION_H
#define IDISTRIBUTION_H

#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>

/**
 * @file IMarginalDistribution.h
 * @brief Interface for one-dimensional marginal probability distributions.
 *
 * This header defines @ref IMarginalDistribution, a minimal runtime-polymorphic
 * interface used by the statistical module to manipulate scalar probability
 * distributions in a uniform way.
 *
 * Implementations are expected to provide:
 * - random sampling,
 * - log-density evaluation,
 * - cumulative distribution function,
 * - inverse cumulative distribution function (quantile / PPF),
 * - first two moments.
 *
 * Typical implementations include:
 * - Gaussian marginals,
 * - flat/uniform marginals,
 * - asymmetric split-Gaussian marginals,
 * - discrete likelihood-based marginals.
 *
 * @see GaussianMarginal
 * @see FlatMarginal
 * @see SplitGaussianMarginal
 * @see LikelihoodMarginal
 */

/**
 * @struct IMarginalDistribution
 * @brief Abstract interface for scalar marginal distributions.
 *
 * A marginal distribution represents a one-dimensional probability law
 * for a nuisance parameter, observable, or latent variable.
 *
 * Implementations should follow these conventions:
 * - @ref rvs returns independent samples,
 * - @ref logpdf returns the natural logarithm of the density (or mass proxy),
 * - @ref cdf returns \f$P(X \le x)\f$,
 * - @ref ppf returns the quantile associated with a probability,
 * - @ref mean and @ref std return the first and second central-moment summary.
 *
 * For distributions where some notions are not naturally defined
 * (for example a discrete likelihood used only as a sampler), implementations
 * may return placeholder values, but this should be clearly documented.
 */
struct IMarginalDistribution {
    /// Virtual destructor for safe polymorphic deletion.
    virtual ~IMarginalDistribution() = default;

    /**
     * @brief Draws random samples from the marginal distribution.
     *
     * @param n Number of samples to generate.
     * @return std::vector<double> of size @p n containing iid draws.
     */
    virtual std::vector<double> rvs(std::size_t n) = 0;

    /**
     * @brief Evaluates the logarithm of the probability density at @p x.
     *
     * For continuous distributions, this is usually \f$\log f(x)\f$.
     * For discrete or sampler-only distributions, the exact semantics depend
     * on the implementation.
     *
     * @param x Point at which the log-density is evaluated.
     * @return Natural logarithm of the density (or implementation-defined proxy).
     */
    virtual double logpdf(double x) = 0;

    /**
     * @brief Evaluates the cumulative distribution function at @p x.
     *
     * @param x Point at which the CDF is evaluated.
     * @return \f$P(X \le x)\f$.
     */
    virtual double cdf(double x) = 0;

    /**
     * @brief Evaluates the quantile function (inverse CDF).
     *
     * @param x Probability in the unit interval.
     * @return Quantile associated with @p x.
     */
    virtual double ppf(double x) = 0;

    /**
     * @brief Returns the mean of the distribution.
     *
     * @return Distribution mean.
     */
    virtual double mean() = 0;

    /**
     * @brief Returns the standard deviation of the distribution.
     *
     * @return Distribution standard deviation.
     */
    virtual double std() = 0;
};

#endif
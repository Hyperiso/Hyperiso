#ifndef COPULATYPE_H
#define COPULATYPE_H

/**
 * @file CopulaType.h
 * @brief Enumerates the supported copula families.
 *
 * This header defines @ref CopulaType, an enum used to select which
 * multivariate copula family should be instantiated by higher-level
 * statistical factories such as @ref CopulaFactory.
 *
 * The copula determines the dependence structure between marginals,
 * independently of the marginal distributions themselves.
 *
 * Currently supported families are:
 *  - Gaussian copula
 *  - Student-t copula
 *
 * @see ICopula
 * @see GaussianCopula
 * @see StudentTCopula
 * @see CopulaFactory
 */

/**
 * @enum CopulaType
 * @brief Identifies the copula family used to model dependence.
 *
 * This enum is typically passed to @ref CopulaFactory::create() together with
 * a compatible configuration object.
 */
enum class CopulaType {
    GAUSSIAN,   ///< Gaussian copula, defined by a correlation matrix.
    STUDENT_T   ///< Student-t copula, defined by a correlation matrix and degrees of freedom.
};

#endif // COPULATYPE_H

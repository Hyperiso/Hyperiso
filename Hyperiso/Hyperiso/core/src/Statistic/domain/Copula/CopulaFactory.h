#ifndef COPULAFACTORY_H
#define COPULAFACTORY_H

#include <memory>
#include <random>
#include <algorithm>
#include <variant>

#include "CopulaType.h"
#include "ICopula.h"
#include "GaussianCopula.h"
#include "StudentTCopula.h"

/**
 * @file CopulaFactory.h
 * @brief Factory utilities for constructing copula objects from typed configurations.
 *
 * This header defines:
 *  - @ref CopulaConfig, a std::variant of supported copula configuration types,
 *  - @ref CopulaFactory, a factory class used to build concrete @ref ICopula
 *    implementations from a @ref CopulaType and a matching configuration.
 *
 * The goal is to decouple:
 *  - the user-facing choice of copula family,
 *  - the concrete implementation class,
 *  - and the configuration payload associated with that family.
 *
 * Typical usage:
 * @code
 *   GaussianCopulaConfig cfg{R};
 *   auto cop = CopulaFactory::create(CopulaType::GAUSSIAN, cfg, 12345);
 *
 *   Vector u = cop->sample_u();
 * @endcode
 *
 * @note The configuration is stored in a std::variant, so only one supported
 *       configuration type can be provided at a time.
 *
 * @see ICopula
 * @see CopulaType
 * @see GaussianCopula
 * @see StudentTCopula
 */

/**
 * @typedef CopulaConfig
 * @brief Variant type gathering all supported copula configuration objects.
 *
 * A @ref CopulaConfig currently holds one of:
 *  - @ref GaussianCopulaConfig
 *  - @ref StudentTCopulaConfig
 *
 * This type is meant to be passed to @ref CopulaFactory::create().
 */
using CopulaConfig = std::variant<GaussianCopulaConfig, StudentTCopulaConfig>;

/**
 * @class CopulaFactory
 * @brief Factory for constructing concrete copulas from a type tag and configuration.
 *
 * CopulaFactory creates heap-allocated concrete implementations of @ref ICopula
 * and returns them through `std::unique_ptr<ICopula>`.
 *
 * Supported constructions currently include:
 *  - @ref GaussianCopula from @ref GaussianCopulaConfig
 *  - @ref StudentTCopula from @ref StudentTCopulaConfig
 *
 * Typical usage:
 * @code
 *   StudentTCopulaConfig cfg;
 *   cfg.R  = corr;
 *   cfg.nu = 5;
 *
 *   auto cop = CopulaFactory::create(CopulaType::STUDENT_T, cfg);
 * @endcode
 *
 * @warning In the current implementation, the actual constructed type is
 *          selected by the active alternative stored in @ref CopulaConfig.
 *          The @p name argument is only used to select a switch branch, but
 *          does not itself enforce consistency with the variant content.
 *          For example, passing `CopulaType::GAUSSIAN` together with a
 *          `StudentTCopulaConfig` still constructs a Student-t copula.
 *
 *          If stricter behavior is desired, the implementation should check
 *          that @p name matches the active type of @p config and throw
 *          otherwise.
 */
class CopulaFactory {
public:
    /**
     * @brief Creates a concrete copula instance.
     *
     * This method builds a concrete @ref ICopula implementation using:
     *  - a user-level copula family tag (@p name),
     *  - a typed configuration payload (@p config),
     *  - and an optional RNG seed (@p seed).
     *
     * Supported pairs are:
     *  - `CopulaType::GAUSSIAN`   with @ref GaussianCopulaConfig
     *  - `CopulaType::STUDENT_T`  with @ref StudentTCopulaConfig
     *
     * @param name   Requested copula family.
     * @param config Configuration object stored in a @ref CopulaConfig variant.
     * @param seed   Optional random seed used to initialize the underlying RNG.
     * @return A unique pointer to the created copula.
     *
     * @throws std::invalid_argument If the copula type is unknown.
     *
     * @warning In the current implementation, no runtime consistency check is
     *          performed between @p name and the active alternative stored in
     *          @p config. The variant content effectively determines the class
     *          that is instantiated.
     */
    static std::unique_ptr<ICopula> create(CopulaType name, CopulaConfig config, unsigned int seed = std::random_device{}());
};

#endif // COPULAFACTORY_H

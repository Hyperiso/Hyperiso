#ifndef WILSON_COEFFICIENT_REGISTRY_H
#define WILSON_COEFFICIENT_REGISTRY_H

#include <functional>
#include <memory>
#include <unordered_map>

#include "GroupDefinition.h"
#include "Wilson.h"

/**
 * @file WilsonCoefficientRegistry.h
 * @brief Factory registry for creating concrete WilsonCoefficient objects by (coefficient, model, backend).
 *
 * This header defines @ref CoefficientRegistry, a lightweight runtime registry used to build
 * the concrete @ref WilsonCoefficient objects that populate a @ref CoefficientGroup.
 *
 * Motivation
 * ----------
 * In the domain layer, a same Wilson coefficient (e.g. C1) can have multiple concrete
 * implementations depending on:
 * - the active physics @ref Model (SM / SUSY / THDM / ...),
 * - the selected computation @ref Backend (Builtin vs Marty).
 *
 * Rather than hard-coding large switch statements in group builders, this registry provides:
 * - explicit registration of "creator" functions (factories),
 * - a single `create()` entry point with well-defined fallback rules.
 *
 * Fallback rules implemented by @ref CoefficientRegistry::create
 * -------------------------------------------------------------
 * When `create(ctx, c)` is called:
 * 1) Try exact match (c, ctx.model, ctx.backend).
 * 2) For a non-SM model using the MARTY backend, route supported SM/BSM
 *    requests through the registered SM/MARTY factory while preserving the
 *    original build context.  The factory then instantiates the target model.
 * 3) Otherwise, if backend is Marty, fallback to Builtin for the same model.
 * 4) If model is not SM, fallback to SM for the same backend.
 *
 * If none match, an exception is thrown.
 *
 * Typical usage
 * -------------
 * @code
 *   CoefficientRegistry reg;
 *   register_B(reg);   // registers creators for the B group coefficients
 *
 *   BuildContext ctx{adapters, Model::SM, Backend::Builtin, ContributionType::SM, gid, "B"};
 *   auto c1 = reg.create(ctx, WCoef::C1);
 * @endcode
 *
 * Registration is typically done per group (see the free functions `register_*`)
 * in a dedicated translation unit.
 *
 * @see BuildContext
 * @see GroupDefinition
 * @see CoefficientGroup
 * @see WilsonCoefficient
 */

/// Shared pointer type used to store coefficients in groups.
using CoefPtr   = std::shared_ptr<WilsonCoefficient>;

/// Factory signature for constructing a coefficient instance given a build context and enum id.
using CoefMaker = std::function<CoefPtr(const BuildContext&, WCoef)>;

/**
 * @class CoefficientRegistry
 * @ingroup DomainModule
 * @brief Registry mapping (WCoef, Model, Backend) -> coefficient factory function.
 *
 * The registry stores creator functions keyed by:
 * - @ref WCoef : which coefficient to create,
 * - @ref Model : which physics model implementation to use,
 * - @ref Backend : built-in computations vs MARTY-generated backend.
 *
 * This is primarily used during the construction of @ref CoefficientGroup objects:
 * each group definition provides membership lists (WCoef) and the registry resolves
 * which concrete @ref WilsonCoefficient subclass to instantiate.
 */
class CoefficientRegistry {
public:
    /**
     * @brief Registers a creator function for a given (coefficient, model, backend).
     *
     * @param c  Coefficient identifier.
     * @param m  Physics model.
     * @param b  Backend choice (Builtin / Marty).
     * @param mk Factory function returning a new coefficient instance.
     */
    void register_creator(WCoef c, Model m, Backend b, CoefMaker mk) {
        table_[key(c,m,b)] = std::move(mk);
    }

    /**
     * @brief Creates a coefficient instance according to the provided build context.
     *
     * The function attempts to find a registered factory using the fallback strategy:
     * - exact match (c, ctx.model, ctx.backend),
     * - non-SM MARTY requests through the SM/MARTY factory with the original context,
     * - Marty -> Builtin fallback for the same model,
     * - model -> SM fallback for the same backend.
     *
     * @param ctx Build context providing model/backend information and adapter access.
     * @param c   Coefficient identifier.
     * @return A newly created coefficient instance.
     * @throws std::runtime_error if no suitable factory is registered.
     */
    CoefPtr create(const BuildContext& ctx, WCoef c) const;
private:
    /**
     * @struct Key
     * @brief Internal key used to index the registry.
     *
     * Uses integer representations of enums to keep hashing simple.
     */
    struct Key {
        int c;
        int m;
        int b;
        bool operator==(const Key& o) const noexcept {
            return c==o.c && m==o.m && b==o.b;
        }
    };

    /**
     * @struct KeyHash
     * @brief Hash functor for @ref Key.
     *
     * Uses a boost-like hash combine to mix the three fields.
     */
    struct KeyHash {
        std::size_t operator()(const Key& k) const noexcept {
            // hash (boost-like)
            std::size_t h = std::hash<int>{}(k.c);
            h ^= std::hash<int>{}(k.m) + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
            h ^= std::hash<int>{}(k.b) + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
            return h;
        }
    };

    /// Builds the internal key from enums.
    static Key key(WCoef c, Model m, Backend b) {
        return Key{(int)c, (int)m, (int)b};
    }

    /// Internal storage of factory functions.
    std::unordered_map<Key, CoefMaker, KeyHash> table_;
};

/**
 * @brief Registers coefficient factories for the B group (C1..C10).
 *
 * Implementations typically register:
 * - Builtin SM/SUSY/THDM concrete classes,
 * - optional Marty creators for SM (and later for other models).
 */
void register_B(CoefficientRegistry&);

/**
 * @brief Registers coefficient factories for the B' group (CP1..CP10, CPQ1, CPQ2).
 */
void register_BPrime(CoefficientRegistry&);

/**
 * @brief Registers coefficient factories for the scalar B group (CQ1, CQ2).
 */
void register_BScalar(CoefficientRegistry&);

/**
 * @brief Registers coefficient factories for charged-current b->c (C_V1_bc, ...).
 */
void register_CC_bc(CoefficientRegistry&);

/**
 * @brief Registers coefficient factories for charged-current b->u (C_V1_bu, ...).
 */
void register_CC_bu(CoefficientRegistry&);

/**
 * @brief Registers coefficient factories for charged-current c->s (C_V1_cs, ...).
 */
void register_CC_cs(CoefficientRegistry&);

/**
 * @brief Registers coefficient factories for charged-current c->d (C_V1_cd, ...).
 */
void register_CC_cd(CoefficientRegistry&);

/**
 * @brief Registers coefficient factories for charged-current s->u (C_V1_su, ...).
 */
void register_CC_su(CoefficientRegistry&);

/**
 * @brief Registers coefficient factories for charged-current d->u (C_V1_du, ...).
 */
void register_CC_du(CoefficientRegistry&);

/**
 * @brief Registers coefficient factories for meson mixing coefficients.
 */
void register_MesonMixing(CoefficientRegistry&);

/**
 * @brief Registers coefficient factories for the K group (CK9, CK10, ...).
 */
void register_K(CoefficientRegistry& reg);

#endif
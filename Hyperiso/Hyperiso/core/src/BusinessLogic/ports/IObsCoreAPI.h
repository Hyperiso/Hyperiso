#ifndef IOBS_CORE_API_H
#define IOBS_CORE_API_H

/**
 * @file IObsCoreAPI.h
 * @brief Minimal observable-layer "core API" concept.
 *
 * This interface is used as a generic adapter to expose small pieces of global
 * runtime state to the observable/business layer (e.g. flags, selected model,
 * feature toggles).
 *
 * It is intentionally minimal: a single getter.
 *
 * @tparam T Type returned by the API (bool, enum, string, ...).
 *
 * @see ObsUseMarty
 */
template <typename T>
class IObsCoreAPI{
public:
    virtual ~IObsCoreAPI() = default;
    /**
     * @brief Returns the current value of the underlying core state.
     */
    inline virtual T get() = 0;
};

#endif
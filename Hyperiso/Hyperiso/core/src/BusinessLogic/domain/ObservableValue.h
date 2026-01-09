#ifndef OBSERVABELEVALUE_H
#define OBSERVABELEVALUE_H

#include <optional>
#include <utility>

#include "Include.h"

/**
 * @struct ObservableValue
 * @brief Container for a computed observable value, optionally binned.
 *
 * This structure represents the numerical result of an observable computation.
 * It stores:
 *
 * - the observable identifier,
 * - the central value,
 * - an optional kinematic bin (e.g. q² range).
 *
 * It is designed to be lightweight and trivially copyable, suitable for:
 * - returning results from observable evaluators,
 * - collecting results in vectors or maps,
 * - passing values to likelihood / fit layers.
 *
 * Typical usage:
 * @code
 * ObservableValue br(BR_B_to_Xs_gamma, 3.36e-4);
 *
 * ObservableValue afb(
 *     AFB_B_to_Kstarmumu,
 *     -0.12,
 *     {1.0, 6.0}   // q^2 bin
 * );
 * @endcode
 *
 * @see ObservableId
 */
struct ObservableValue {
    /// Identifier of the observable.
    ObservableId id;

    /// Central numerical value of the observable.
    double value;

    /// Optional bin definition (lower, upper).
    std::optional<std::pair<double, double>> bin;

    /**
     * @brief Construct an unbinned observable value.
     * @param id    Observable identifier.
     * @param value Central value.
     */
    ObservableValue(ObservableId id, double value) : id(id), value(value) {};

    /**
     * @brief Construct a binned observable value.
     * @param id    Observable identifier.
     * @param value Central value.
     * @param bin   Kinematic bin (lower, upper).
     */
    ObservableValue(ObservableId id, double value, std::pair<double, double> bin) : ObservableValue(id, value) { this->bin = bin; };
};

#endif // OBSERVABELEVALUE_H

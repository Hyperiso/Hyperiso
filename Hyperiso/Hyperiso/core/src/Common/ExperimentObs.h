#ifndef EXPERIMENT_OBS_H
#define EXPERIMENT_OBS_H

#include "Include.h"

/**
 * @file ExperimentObs.h
 * @brief Experiment-scoped observable identifiers.
 *
 * This header defines @ref ExperimentObs, a small value type combining an
 * experiment label with a @ref BinnedObservableId. It is intended to uniquely
 * identify a measured observable in the context of a specific experiment.
 */

/**
 * @struct ExperimentObs
 * @brief Identifies an observable in the context of a specific experiment.
 *
 * ExperimentObs combines the name of an experiment with a binned observable
 * identifier. It can distinguish measurements of the same observable performed
 * by different experiments, or measurements belonging to different
 * experimental contexts.
 *
 * The observable part is stored as a @ref BinnedObservableId. Convenience
 * constructors are provided for already-binned observables, plain
 * @ref ObservableId values, and @ref Observables enum values. In the latter
 * two cases, the bin boundaries are initialized to {0.0, 0.0}.
 *
 * The type supports equality comparison, lexicographic ordering, hashing, and
 * human-readable string formatting.
 */
struct ExperimentObs {
    std::string experiment; /**< Name or label of the experiment. */
    BinnedObservableId obs; /**< Observable identifier, including optional bin information. */

    /**
     * @brief Default constructor.
     *
     * Leaves the experiment name and observable identifier default-initialized.
     */
    ExperimentObs() = default;

    /**
     * @brief Constructs an experiment-specific binned observable.
     *
     * @param exp Name or label of the experiment.
     * @param o Binned observable identifier.
     */
    ExperimentObs(std::string exp, BinnedObservableId o)
        : experiment(std::move(exp)), obs(std::move(o)) {}

    /**
     * @brief Constructs an experiment-specific unbinned observable.
     *
     * The observable is wrapped into a BinnedObservableId with bin boundaries
     * set to {0.0, 0.0}.
     *
     * @param exp Name or label of the experiment.
     * @param o Observable identifier.
     */
    ExperimentObs(std::string exp, ObservableId o)
        : experiment(std::move(exp)), obs(BinnedObservableId(o, {0.,0.})) {}

    /**
     * @brief Constructs an experiment-specific observable from an enum value.
     *
     * The observable enum value is converted to an ObservableId through
     * @ref ObservableMapper::to_id and wrapped into a BinnedObservableId with
     * bin boundaries set to {0.0, 0.0}.
     *
     * @param exp Name or label of the experiment.
     * @param o Observable enum value.
     */
    ExperimentObs(std::string exp, Observables o)
        : experiment(std::move(exp)), obs(BinnedObservableId(ObservableMapper::to_id(o), {0.,0.})) {}

    /**
     * @brief Equality comparison operator.
     *
     * Two ExperimentObs objects are equal if and only if both their experiment
     * labels and their binned observable identifiers are equal.
     *
     * @param other Object to compare with.
     * @return true if both objects identify the same observable in the same
     *         experiment, false otherwise.
     */
    bool operator==(ExperimentObs const& other) const noexcept {
        return experiment == other.experiment
            && obs == other.obs;
    }

    /**
     * @brief Strict weak ordering for experiment-specific observables.
     *
     * Ordering is lexicographic on the experiment label and then on the binned
     * observable identifier. This allows the type to be used as a key in
     * ordered containers such as std::map or std::set.
     *
     * @param other Object to compare with.
     * @return true if this object is ordered before @p other.
     */
    bool operator<(ExperimentObs const& other) const noexcept {
        return std::tie(experiment, obs)
             < std::tie(other.experiment, other.obs);
    }

    /**
     * @brief Returns a human-readable representation of the object.
     *
     * @return String of the form "experiment :: observable [low, high]".
     */
    std::string str() const {
        std::stringstream ss;
        ss << experiment << " :: " << obs.str();
        return ss.str();
    }
};

/**
 * @brief Stream output operator for ExperimentObs.
 *
 * Writes the human-readable representation of the experiment-specific
 * observable to the output stream.
 *
 * @param os Output stream.
 * @param x Experiment-specific observable to print.
 * @return Reference to the output stream.
 */
inline std::ostream& operator<<(std::ostream& os, ExperimentObs const& x) {
    os << x.experiment << " :: " << x.obs.str();
    return os;
}

/**
 * @brief Hash functor specialization for ExperimentObs.
 *
 * This specialization allows ExperimentObs to be used as a key in
 * std::unordered_map or std::unordered_set. The hash combines the experiment
 * label and the hash of the binned observable identifier.
 */
template<>
struct std::hash<ExperimentObs> {
    /**
     * @brief Computes the hash of an experiment-specific observable.
     *
     * @param x ExperimentObs object to hash.
     * @return Hash value for @p x.
     */
    std::size_t operator()(ExperimentObs const& x) const noexcept {
        std::size_t h = std::hash<std::string>{}(x.experiment);

        auto mix = [](std::size_t& seed, std::size_t v) {
            seed ^= v + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        };

        mix(h, std::hash<BinnedObservableId>{}(x.obs));
        return h;
    }
};

#endif // EXPERIMENT_OBS_H

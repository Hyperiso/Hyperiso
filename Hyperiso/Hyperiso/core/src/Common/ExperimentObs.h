#ifndef EXPERIMENT_OBS_H
#define EXPERIMENT_OBS_H

#include "Include.h"

struct ExperimentObs {
    std::string experiment;
    BinnedObservableId obs;

    ExperimentObs() = default;

    ExperimentObs(std::string exp, BinnedObservableId o)
        : experiment(std::move(exp)), obs(std::move(o)) {}

    ExperimentObs(std::string exp, ObservableId o)
        : experiment(std::move(exp)), obs(BinnedObservableId(o, {0.,0.})) {}

    ExperimentObs(std::string exp, Observables o)
        : experiment(std::move(exp)), obs(BinnedObservableId(ObservableMapper::to_id(o), {0.,0.})) {}

    bool operator==(ExperimentObs const& other) const noexcept {
        return experiment == other.experiment
            && obs == other.obs;
    }

    bool operator<(ExperimentObs const& other) const noexcept {
        return std::tie(experiment, obs)
             < std::tie(other.experiment, other.obs);
    }

    std::string str() const {
        std::stringstream ss;
        ss << experiment << " :: " << obs.str();
        return ss.str();
    }
};

inline std::ostream& operator<<(std::ostream& os, ExperimentObs const& x) {
    os << x.experiment << " :: " << x.obs.str();
    return os;
}

template<>
struct std::hash<ExperimentObs> {
    std::size_t operator()(ExperimentObs const& x) const noexcept {
        std::size_t h = std::hash<std::string>{}(x.experiment);

        auto mix = [](std::size_t& seed, std::size_t v) {
            seed ^= v + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        };

        mix(h, std::hash<BinnedObservableId>{}(x.obs));
        return h;
    }
};

#endif
#include "Include.h"

struct ExperimentObs {
    std::string experiment;
    BinnedObservableId obs;

    ExperimentObs() = default;

    ExperimentObs(std::string exp, BinnedObservableId o)
        : experiment(std::move(exp)), obs(std::move(o)) {}

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
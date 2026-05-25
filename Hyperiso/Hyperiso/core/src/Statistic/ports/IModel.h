#ifndef IMODEL_H
#define IMODEL_H

#include <cstddef>
#include <vector>
#include <map>
#include "Include.h"
#include "ObservableValue.h"

using Vec = std::vector<double>;

class IModel {
public:
    virtual ~IModel() = default;

    virtual std::map<ObservableId, std::vector<ObservableValue>> predict_optimized(const std::map<ParamId, double>& p, const std::map<ParamId, double>& eta) = 0;
    virtual std::size_t n_observables() const = 0;
    virtual std::unordered_set<ParamId> get_obs_deps(ObservableId id) = 0;
    virtual std::vector<BinnedObservableId> get_obs_ids() = 0;
    virtual void compute_observables() const = 0;
};

#endif
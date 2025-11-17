#ifndef IMODEL_H
#define IMODEL_H

#include <cstddef>
#include <vector>


using Vec = std::vector<double>;


class IModel {
public:
virtual ~IModel() = default;

virtual Vec predict(const Vec& p, const Vec& eta) = 0;
virtual std::size_t n_observables() const = 0;
virtual void add_observables(std::map<ObservableId, QCDOrder> obs_ids) = 0;
virtual std::unordered_set<ParamId> get_obs_deps(ObservableId id) = 0;
};

#endif
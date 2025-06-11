#ifndef COMPOUND_H
#define COMPOUND_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include "ObsParameterProxy.h"
#include "ObsParameterMutator.h"
#include "CorrelationProxy.h"

struct Estimate {
    scalar_t central_value;
    scalar_t stat_std;
    scalar_t syst_std;

    scalar_t combined_std() {
        return std::hypot(stat_std, syst_std);
    }

    friend std::ostream& operator<<(std::ostream& os, Estimate est) {
        os << est.central_value << " +- " << est.stat_std << " (stat)" << " +- " << est.syst_std << "(syst)";
        return os;
    }
};

class Compound {

protected:
    std::unordered_set<ParamId> dependences;
    std::unordered_map<ParamId, scalar_t> gradient;
    scalar_t central_value {NAN};

    scalar_t compute_pdv(const ParamId& param_name) const;
    std::unordered_set<ParamId> get_common_dependences_with(const Compound& other) const;

public:
    virtual scalar_t eval() const = 0;
    void add_dependence(const ParamId& param_name);
    void add_dependences(const std::unordered_set<ParamId>& param_names);
    void update_gradient();
    const std::unordered_set<ParamId>& get_dependences() const;
    const std::unordered_map<ParamId, scalar_t>& get_gradient() const;
    scalar_t variance();
    scalar_t correlation_with(const Compound& other) const;
    const std::unordered_map<ParamId, scalar_t> get_leading_uncertainties(size_t n) const;
    const std::unordered_map<ParamId, scalar_t> get_uncertainties() const;
    Estimate get_estimate();
    void print_gradient(std::ostream& os) const;

};

#endif // __COMPOUND_H__
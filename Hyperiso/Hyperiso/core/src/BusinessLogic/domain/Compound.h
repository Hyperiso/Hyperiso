#ifndef COMPOUND_H
#define COMPOUND_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include "Parameter.h"
// #include "CorrelationRepo.h"
#include "ObsParameterProxy.h"
#include "ObsParameterMutator.h"
#include "CorrelationProxy.h"

class Compound {

protected:
    std::vector<ParamId> dependences;
    std::map<ParamId, double> gradient;
    double central_value {NAN};

    double compute_pdv(const ParamId& param_name) const;
    std::vector<ParamId> get_common_dependences_with(const Compound& other) const;

public:
    virtual double eval() const = 0;
    void add_dependence(const ParamId& param_name);
    void add_dependences(const std::vector<ParamId>& param_names);
    void update_gradient();
    const std::vector<ParamId>& get_dependences() const;
    const std::map<ParamId, double>& get_gradient() const;
    double variance();
    double correlation_with(const Compound& other) const;
    const std::map<ParamId, double> get_leading_uncertainties(size_t n) const;
    const std::map<ParamId, double> get_uncertainties() const;
    void print_gradient(std::ostream& os) const;

};

#endif // __COMPOUND_H__
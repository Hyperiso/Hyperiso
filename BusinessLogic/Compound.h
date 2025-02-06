#ifndef __COMPOUND_H__
#define __COMPOUND_H__

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "Parameter.h"
#include "Matrix.h"

class Compound {

protected:
    std::vector<ParamId> dependences;
    SparseMatrix<ParamId> param_corr;
    std::map<ParamId, double> gradient;

    void read_param_covariance();
    double compute_pdv(const ParamId& param_name) const;
    void update_gradient();
    std::vector<ParamId> get_common_dependences_with(const Compound& other) const;

public:
    virtual double eval() const = 0;
    void add_dependence(const ParamId& param_name);
    void add_dependences(const std::vector<ParamId>& param_names);
    const std::vector<ParamId>& get_dependences() const;
    const std::map<ParamId, double>& get_gradient() const;
    double variance() const;
    double correlation_with(const Compound& other) const;
    const std::map<ParamId, double> get_leading_uncertainties(size_t n) const;
    const std::map<ParamId, double> get_uncertainties() const;
    void print_gradient(std::ostream& os) const;

};

#endif // __COMPOUND_H__
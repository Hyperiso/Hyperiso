#ifndef MODELEVALUATOR_H
#define MODELEVALUATOR_H

#include "Compound.h"
#include "Observable.h"
#include "General.h"
#include "Matrix.h"
#include "CorrelationProxy.h"
#include <vector>
#include <string>

class ModelEvaluator {

private:

    std::vector<std::shared_ptr<Observable>> observables;
    SparseMatrix<Observables> th_cov_mtx;
    SparseMatrix<Observables> exp_cov_mtx;

    void update_th_covariance();
    // void update_exp_data(); //TODO: correlation proxy

public:

    ModelEvaluator();
    ModelEvaluator(const std::vector<std::shared_ptr<Observable>>& observables);

    bool has_observable(Observables id);
    void add_observable(std::shared_ptr<Observable> obs);
    void remove_observable(Observables id);
    std::shared_ptr<Observable> find_from_id(Observables id);
    SparseMatrix<Observables> get_covariance();
    double chi2();

};

#endif // __MODELEVALUATOR_H__

#ifndef MODELEVALUATOR_H
#define MODELEVALUATOR_H

#include "Compound.h"
#include "Observable.h"
#include "Include.h"
#include "Matrix.h"
#include "CorrelationProxy.h"
#include <vector>
#include <string>

class ModelEvaluator {

private:

    std::unordered_map<ObservableId, std::shared_ptr<Observable>> observables;
    SparseMatrix<ObservableId> th_cov_mtx;
    SparseMatrix<ObservableId> exp_cov_mtx;

    void update_th_covariance();
    void update_exp_covariance();

public:

    ModelEvaluator();
    ModelEvaluator(const std::unordered_set<std::shared_ptr<Observable>>& observables);

    bool has_observable(Observables id);
    bool has_observable(ObservableId id);

    void add_observable(std::shared_ptr<Observable> obs);
    void remove_observable(Observables id);
    void remove_observable(ObservableId id);

    SparseMatrix<ObservableId> get_covariance();
    double chi2();

};

#endif // __MODELEVALUATOR_H__

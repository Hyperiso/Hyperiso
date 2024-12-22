#ifndef __MODELEVALUATOR_H__
#define __MODELEVALUATOR_H__

#include "Compound.h"
#include "Observable.h"
#include "General.h"
#include "Matrix.h"
#include <vector>
#include <string>

class ModelEvaluator {

private:

    std::vector<std::shared_ptr<Observable>> observables;
    SparseMatrix<Observables> th_cov_mtx;
    SparseMatrix<Observables> exp_cov_mtx;

    void update_th_covariance();
    void update_exp_data();

public:

    ModelEvaluator();
    ModelEvaluator(const std::vector<std::shared_ptr<Observable>>& observables);

    std::shared_ptr<Observable> find_from_id(Observables id);
    SparseMatrix<Observables> get_covariance() const;
    double chi2() const;

};

#endif // __MODELEVALUATOR_H__

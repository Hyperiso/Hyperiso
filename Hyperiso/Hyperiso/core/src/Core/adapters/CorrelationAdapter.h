#ifndef __CORRELATIONCREATOR_H__
#define __CORRELATIONCREATOR_H__

#include "General.h"
#include "Math.h"
#include "DBNode.h"
#include "ParameterRouter.h"
#include "CorrelationRepo.h"

class CorrelationAdapter {
public:
    template<typename T>
    static CorrelationMatrixPair<T> from_db_node(std::shared_ptr<Node> root);

private:
    static void emplace_correlation(CorrelationMatrixPair<ParamId>& corr_matrices, std::shared_ptr<Node> leaf);
    static void emplace_correlation(CorrelationMatrixPair<Observables>& corr_matrices, std::shared_ptr<Node> leaf);
};

#endif // __CORRELATIONCREATOR_H__

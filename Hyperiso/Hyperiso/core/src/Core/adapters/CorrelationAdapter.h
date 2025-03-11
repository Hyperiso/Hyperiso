#ifndef __CORRELATIONCREATOR_H__
#define __CORRELATIONCREATOR_H__

#include "General.h"
#include "ParameterRouter.h"
#include "IDataLoader.h"

template<typename T>
class CorrelationLoader : public IDataLoader<CorrelationMatrixPair<T>> {
public:
    void load(std::shared_ptr<CorrelationMatrixPair<T>> dest, fs::path src_file) override;

private:
    static void emplace_correlation(std::shared_ptr<CorrelationMatrixPair<ParamId>> corr_matrices, std::shared_ptr<Node> leaf);
    static void emplace_correlation(std::shared_ptr<CorrelationMatrixPair<Observables>> corr_matrices, std::shared_ptr<Node> leaf);
};

#endif // __CORRELATIONCREATOR_H__

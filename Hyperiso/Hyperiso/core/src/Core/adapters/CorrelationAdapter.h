#ifndef CORRELATIONCREATOR_H
#define CORRELATIONCREATOR_H

#include "Include.h"
#include "ParameterRouter.h"
#include "IDataLoader.h"
#include "CorrelationRepo.h"
#include "DBNode.h"
#include "DBNodeProviderFactory.h"

/**
 * @file CorrelationLoader.h
 * @brief Loads correlation matrices from a file into CorrelationMatrixPair structures.
 *
 * This file defines the CorrelationLoader class, which specializes the IDataLoader
 * interface to load statistical and systematic correlation matrices.
 */

/**
 * @class CorrelationLoader
 * @ingroup DataLoadersModule
 * @brief Loader class to populate CorrelationMatrixPair from a file.
 *
 * @tparam T Type of keys for the correlation matrix (e.g., ParamId or ObservableId).
 */
template<typename T>
class CorrelationLoader : public IDataLoader<CorrelationMatrixPair<T>> {
public:

    /**
     * @brief Loads correlation matrices from a source file into a destination CorrelationMatrixPair.
     *
     * @param dest     Shared pointer to destination CorrelationMatrixPair.
     * @param src_file Path to the source file.
     */
    void load(std::shared_ptr<CorrelationMatrixPair<T>> dest, fs::path src_file) override;

private:
    /**
     * @brief Inserts a correlation entry into a correlation matrix.
     *
     * @param corr_matrices Shared pointer to the correlation matrix pair.
     * @param leaf          DBNode containing correlation data.
     */
    static void emplace_correlation(std::shared_ptr<CorrelationMatrixPair<T>> corr_matrices, std::shared_ptr<DBNode> leaf);
};

#endif // CORRELATIONCREATOR_H

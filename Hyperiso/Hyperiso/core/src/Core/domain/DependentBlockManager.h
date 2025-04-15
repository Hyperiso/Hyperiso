#ifndef DEPENDENT_BLOCK_MANAGER_H
#define DEPENDENT_BLOCK_MANAGER_H

#include "Block.h"
#include <map>
#include <unordered_map>
#include <memory>
#include <vector>
#include <iostream>

class DependentBlockManager {
public:
    /**
     * @brief Ajoute un DependentBlock avec plusieurs sources provenant de différentes instances de Parameters.
     * @param name Nom du DependentBlock.
     * @param sources Map des sources (clé = nom du block, valeur = shared_ptr vers le block).
     * @param recalculateFunc Fonction de recalcul du block dépendant.
     */
    static void addDependentBlock(
        const std::string& name,
        const std::unordered_map<ParameterType, std::vector<std::string>>& source_names,
        ParameterType dest,
        std::function<void(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>)> recalculateFunc
    );

    /**
     * @brief Ajoute un DependentParameter avec plusieurs sources.
     * @param name Id du DependentBlock.
     * @param sources Set des sources.
     * @param recalculateFunc Fonction de recalcul du paramètre.
     */
    static void addDependentParameter(
        const ParamId &pid,
        const std::unordered_set<ParamId> &source_pids,
        DepParamUpdateFunc recalculateFunc
    );

    static void removeDependentBlock(const std::string& name, ParameterType src);

    static void update(const std::string& name, ParameterType src);

private:
    DependentBlockManager() = default;
    ~DependentBlockManager() = default;
    DependentBlockManager(const DependentBlockManager&) = delete;
    DependentBlockManager& operator=(const DependentBlockManager&) = delete;
};

#endif // DEPENDENT_BLOCK_MANAGER_H
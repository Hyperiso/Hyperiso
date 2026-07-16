#ifndef ILHA_PROTOTYPE_REGISTRY_H
#define ILHA_PROTOTYPE_REGISTRY_H

#include <cstddef>
#include <vector>

#include "Include.h"

/**
 * @brief Value object describing an additional LHA block prototype.
 */
struct LhaPrototypeSpec {
    BlockName blockName;
    size_t itemCount = 2;
    size_t valueIdx = 1;
    int scaleIdx = -1;
    int rgIdx = -1;
    int binIdx = -1;
    bool globalScale = false;
};

/**
 * @brief Port used by the domain to register LHA parser prototypes.
 *
 * Implementations live outside MemoryManager. MemoryManager only depends on
 * this abstraction, while HyperisoMaster decides which concrete adapter is
 * injected.
 */
class ILhaPrototypeRegistry {
public:
    virtual ~ILhaPrototypeRegistry() = default;

    virtual void add_lha_prototype(BlockName blockName,
                                   size_t itemCount = 2,
                                   size_t valueIdx = 1,
                                   int scaleIdx = -1,
                                   int rgIdx = -1,
                                   int binIdx = -1,
                                   bool globalScale = false) = 0;

    virtual void add_lha_prototypes(const std::vector<LhaPrototypeSpec>& prototypes)
    {
        for (const auto& prototype : prototypes) {
            add_lha_prototype(prototype.blockName,
                              prototype.itemCount,
                              prototype.valueIdx,
                              prototype.scaleIdx,
                              prototype.rgIdx,
                              prototype.binIdx,
                              prototype.globalScale);
        }
    }
};

#endif // ILHA_PROTOTYPE_REGISTRY_H

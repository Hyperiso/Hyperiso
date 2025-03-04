#ifndef __IDATALOADER_H__
#define __IDATALOADER_H__

#include <memory>
#include "BlockAccessor.h"
#include "CorrelationRepo.h"
#include <any>
#include <typeindex>

class IDataLoader {
public:
    virtual std::shared_ptr<BlockAccessor> load_param_blocks() = 0;
    virtual std::any load_correlations(std::type_index type) = 0;

    virtual ~IDataLoader() = default;
};

#endif // __IDATALOADER_H__

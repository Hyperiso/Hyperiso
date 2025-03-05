#ifndef __IDATALOADER_H__
#define __IDATALOADER_H__

#include <memory>
#include "BlockAccessor.h"
#include "CorrelationRepo.h"
#include <any>
#include <typeindex>

template<typename T>
class IDataLoader {
public:
    virtual void load(T* dest, shared_ptr<Node> src) = 0;

    virtual ~IDataLoader() = default;
};

#endif // __IDATALOADER_H__

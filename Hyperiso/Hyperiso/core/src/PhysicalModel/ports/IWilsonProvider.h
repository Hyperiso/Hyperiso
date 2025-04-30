#ifndef __IWILSONPROVIDER_H__
#define __IWILSONPROVIDER_H__

#include "Math.h"
#include "AbstractConfig.h"
#include <memory>

class IWilsonProvider {
public:
    virtual ~IWilsonProvider() = default;

    virtual scalar_t get(std::shared_ptr<AbstractConfig>) = 0;
};

#endif // __IWILSONPROVIDER_H__

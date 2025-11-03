#ifndef IWILSONPROVIDER_H
#define IWILSONPROVIDER_H

#include "Math.h"
#include "AbstractConfig.h"
#include <memory>
#include "IWilsonBuilder.h"
#include "Include.h"

template<typename BuilderType>
class IWilsonProvider {
public:
    virtual ~IWilsonProvider() = default;

    virtual scalar_t get(std::shared_ptr<AbstractConfig>) = 0;
    virtual std::unordered_set<WilsonBasis> get_bases(WGroupId) = 0;
    virtual std::shared_ptr<BuilderType> get_builder() = 0;
};

#endif // IWILSONPROVIDER_H

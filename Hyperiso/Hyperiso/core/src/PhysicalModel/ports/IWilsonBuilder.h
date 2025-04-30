#ifndef __IWILSONBUILDER_H__
#define __IWILSONBUILDER_H__

#include <memory>
#include "IWilsonProvider.h"

template<typename ConfigType, typename GroupIdType, typename ProviderType>
class IWilsonBuilder {
public:
    virtual ~IWilsonBuilder() = default;

    virtual void build(ConfigType) = 0;
    virtual void add(ConfigType) = 0;
    virtual void switch_basis(GroupIdType) = 0;
    virtual std::shared_ptr<ProviderType> get_wilson_provider() = 0;
};

#endif // __IWILSONBUILDER_H__

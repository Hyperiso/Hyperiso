#ifndef IWILSONBUILDER_H
#define IWILSONBUILDER_H

#include <memory>
#include "IWilsonProvider.h"

template<typename ConfigType, typename ProviderType>
class IWilsonBuilder {
public:
    virtual ~IWilsonBuilder() = default;

    virtual void build(ConfigType) = 0;
    virtual void add(ConfigType) = 0;
    virtual std::shared_ptr<ProviderType> get_wilson_provider() = 0;
};

#endif // IWILSONBUILDER_H

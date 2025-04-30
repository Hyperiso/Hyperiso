#ifndef __IOBSWILSONBUILDER_H__
#define __IOBSWILSONBUILDER_H__

#include "AbstractConfig.h"
#include <memory>

template<typename ProxyType, typename IdType>
class IObsWilsonBuilder {
public: 
    virtual ~IObsWilsonBuilder() = default;	

    virtual void build(std::shared_ptr<AbstractConfig>) = 0;
    virtual void switch_basis(IdType) = 0;
    virtual std::shared_ptr<ProxyType> get_proxy() = 0;
};

#endif // __IOBSWILSONBUILDER_H__

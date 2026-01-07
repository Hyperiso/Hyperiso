#ifndef IOBSWILSONBUILDER_H
#define IOBSWILSONBUILDER_H

#include "AbstractConfig.h"
#include <memory>

class IObsWilsonProxy;

class IObsWilsonBuilder {
public: 
    virtual ~IObsWilsonBuilder() = default;	

    virtual void build(std::shared_ptr<AbstractConfig>) = 0;

    virtual std::shared_ptr<IObsWilsonProxy> get_proxy() = 0;
};

#endif // IOBSWILSONBUILDER_H

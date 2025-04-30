#ifndef __IOBSWILSONBUILDER_H__
#define __IOBSWILSONBUILDER_H__

template<typename InterfaceType, typename ConfigType>
class IObsWilsonBuilder {
public: 
    virtual ~IObsWilsonBuilder() = default;	

    virtual void build(ConfigType) = 0;
    virtual std::shared_ptr<InterfaceType> interface() = 0;
};

#endif // __IOBSWILSONBUILDER_H__

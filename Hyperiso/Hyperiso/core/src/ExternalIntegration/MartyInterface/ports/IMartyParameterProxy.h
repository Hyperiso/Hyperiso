#ifndef __IMARTYPARAMETERPROXY_H__
#define __IMARTYPARAMETERPROXY_H__

#include "scalar.h"

template<typename T, typename V>
class IMartyParameterProxy {
public:
    virtual scalar_t operator()(const T& x, const V& y) const = 0 ;
};

#endif

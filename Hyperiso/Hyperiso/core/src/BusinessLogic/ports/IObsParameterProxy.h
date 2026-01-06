#ifndef IOBS_PARAMETER_PROXY_H
#define IOBS_PARAMETER_PROXY_H

#include <concepts>

#include "Math.h"
#include "ParameterProvider.h"

template<typename T, typename V, typename U, typename W>
class IObsParameterProxy {
public:
    virtual scalar_t operator()(const T& pid, V d_type) = 0;

    virtual scalar_t operator()(const U& block, const W& id, V d_type) const = 0;

    virtual std::shared_ptr<Parameter> get_parameter(const T& pid) const = 0;
};

#endif